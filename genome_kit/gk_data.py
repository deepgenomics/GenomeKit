# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
""" This module automatically fetches GenomeKit data.

A user should first use :py:func:`~genome_kit.gk_data.upload_file` to upload a
file to make it available other GenomeKit users.

After the upload, the file can be used by any GenomeKit function via
:py:func:`~genome_kit.gk_data.get_file`, which will download the file on demand
and return its local path.

Example
-------
>>> upload_file('/local/path/hg38.2bit', 'hg38.2bit')
>>> get_file('hg38.2bit')
"/Users/example/Application Support/genome_kit/hg38.2bit"

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import logging
import os
from typing import Dict, List
from importlib_metadata import EntryPoint, entry_points


from ._util import makedirs
from ._gk_data_config import _config
from . import _cxx
from .data_manager import DataManager, DefaultDataManager, ProgressPercentage
import tempfile
import ftplib
import hashlib
import calendar
import time

from urllib.parse import urlparse
from urllib.request import urlopen


eps = entry_points(group="genomekit.plugins.data_manager")
list_eps: List[EntryPoint] = list(eps)
len_eps = len(list_eps)
try:
    if len_eps > 0:
        DataManagerImpl = list_eps[0].load()
    if len_eps > 1:
        logging.info("Multiple data manager plugins found. "
                     f"Using the first one: {list_eps[0].name}")
    if len_eps == 0:
        DataManagerImpl = DefaultDataManager
except Exception as e:
    logging.debug(e, exc_info=True)
    if len(eps) > 0:
        logging.warning("Failed to load the data manager plugin. "
                        "Falling back to the default data manager.")
    DataManagerImpl = DefaultDataManager
# Register an implementation of DataManager if you need to upload
# new data
data_manager: DataManager = DataManagerImpl(_config["DATA_DIR"])


def get_file(filename):
    """compatibility wrapper for :py:meth:`~genome_kit.DataManager.get_file`"""
    return data_manager.get_file(filename)


def upload_file(filepath:str, filename:str, metadata:Dict[str, str]=None):  # pragma: no cover
    """compatibility wrapper for :py:meth:`~genome_kit.DataManager.upload_file`"""

    chrom_sizes_ext = ".chrom.sizes"
    if filename.endswith(chrom_sizes_ext):
        refg_name = filename[:-len(chrom_sizes_ext)]
        hashval = _cxx.Genome._refg_hash(refg_name)
        print(f"Detected upload of a refg. Creating and uploading a hash lookup file {hashval}.hash")
        tmpfilename = None
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write(refg_name)
            tmpfilename = f.name
        data_manager.upload_file(tmpfilename, f"{hashval}.hash", metadata)

    return data_manager.upload_file(filepath, filename, metadata)

@_cxx.register
def resolve_datafile_path(path):
    """Convert an abstract file path into a concrete, possibly resolving to a GenomeKit data file.

    If path is a regular path, then it is simply returned.

    If path contains ```'{GENOMEKIT_DATA_DIR}'```, then the file is presumed to reside in the GenomeKit
    data file repository. In that case, `get_file` will be attempted on the file name, so that
    the file may be downloaded and versioned if necessary. The path to the final file is returned.
    """
    # If this path includes the special GENOMEKIT_DATA_DIR signifier,
    # process the path and try to download the file if missing.
    # Otherwise just return the file name as-is.
    if "{GENOMEKIT_DATA_DIR}" in path:  # pragma: no cover
        rel_path = path.replace("{GENOMEKIT_DATA_DIR}", "").lstrip('/\\')
        path = get_file(rel_path)
    return path


def _try_ftp_last_modified_date(url):
    """Tries to determine the last-modified date, as UTC in seconds since epoch,
    of a file using the FTP LIST command. If not date could be determined
    using the LIST method, returns None.
    Raises an error if the url does not exist or if it is not a file.
    """

    # Parse URL and check that it's FTP
    url = urlparse(url)
    if url.scheme != 'ftp':
        return None

    # Open an FTP connection
    ftp = ftplib.FTP(url.netloc, user='anonymous', timeout=20)

    # Run the MDTM command to get time in "YYYYMMDDhhmmss" format where MM is 01-12, DD is 01-31, hh is 00-23
    last_modified_str = ftp.sendcmd("MDTM " + url.path).split()[-1]
    try:
        ftp.quit()
    except (OSError, EOFError):
        pass
    finally:
        if ftp.sock is not None:
            ftp.close()

    # Convert to UTC as epoch
    timestruct = time.strptime(last_modified_str, "%Y%m%d%H%M%S")
    return calendar.timegm(timestruct)


def wget(url, dst=None, timestamping=False, progress=False):
    """Download a file from the URL.

    If `dst` is unspecified, the file will be downloaded to
    the system temp directory and given the same name.

    If the last-modified time of the remote file is available, it
    will be assigned to the downloaded file.

    Parameters
    ----------
    url : :py:class:`str`
        The URL of the source file.
    dst : :py:class:`str`
        The path to the local destination.
    timestamping: :py:class:`bool`
        If true, only download when last-modified date if the remote file
        is newer than the local file. If the server does not provide the
        last-modified date, the local file will be kept by default.
    progress: :py:class:`bool`
        If true, show a progress bar.

    Returns
    -------
    :py:class:`str`
        The path to the local version of the file.

    Raises
    ------
    :py:exc:`IOError`
        There was an error opening `dst`.
    :py:exc:`HTTPError`
        There was an error opening the URL.
    :py:exc:`socket.timeout`
        The connection timed out or was interrupted.
    """

    if not isinstance(url, str):
        raise TypeError("url must be str")

    # Get path to output file.
    if dst is None:
        # Use an auto-generated filename based on the URL filename plus a hash of full URL
        tempdir = tempfile.gettempdir()
        filename = "genomekit.{}.{}".format(hashlib.sha1(url.encode('ascii')).hexdigest(), url.split('/')[-1])
        dst = os.path.join(tempdir, filename)

    elif not isinstance(dst, str):
        raise TypeError("dst must be str")

    # Name of partial download file
    dst_part = dst + ".part"

    # Open the URL and check header for timestamp
    with urlopen(url, timeout=20) as request:
        # Fetch remote last-modified time
        remote_last_modified = request.headers['last-modified']
        if remote_last_modified:  # UTC time
            remote_last_modified = calendar.timegm(time.strptime(remote_last_modified, '%a, %d %b %Y %H:%M:%S %Z'))
        else:
            remote_last_modified = _try_ftp_last_modified_date(url)

        # Fetch remote file size
        remote_file_size = request.headers.get("Content-Length", None)
        if remote_file_size is not None:
            remote_file_size = int(remote_file_size)

        # If dst exists, only proceed with download if can determine that remote file is newer
        if timestamping and os.path.exists(dst):
            local_last_modified = os.path.getmtime(dst)
            if remote_last_modified is None or local_last_modified >= remote_last_modified:
                return dst

        # Download the file content
        try:
            makedirs(os.path.dirname(dst))
            with open(dst_part, 'wb') as f:

                # Create a progress bar requested
                if progress:
                    progress = ProgressPercentage(dst, remote_file_size, '[GenomeKit] downloading')  # pragma: no cover

                while True:
                    # Read a chunk into dst
                    buffer = request.read(2**18)
                    if not buffer:
                        break
                    f.write(buffer)

                    # Report progress
                    if progress:
                        progress(len(buffer))  # pragma: no cover

        finally:
            # Set the modification time on the file. This should happen whether
            # the download completed successfully or not.
            # Skip codecov because don't currently know of URL that triggers this case.
            if remote_last_modified is not None:  # pragma: no cover
                os.utime(dst_part, (remote_last_modified, remote_last_modified))

    if os.path.isfile(dst):
        os.unlink(dst)
    os.rename(dst_part, dst)

    return dst
