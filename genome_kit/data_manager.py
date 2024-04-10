# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.

import base64
import hashlib
import logging
import os
import tempfile
from abc import ABC
from contextlib import contextmanager
from functools import wraps, lru_cache
from pathlib import Path
from typing import Dict

from google.cloud import storage
from tqdm.auto import tqdm
from tqdm.utils import ObjectWrapper

logger = logging.getLogger(__name__)

class ProgressPercentage(object):  # pragma: no cover
    """
    Display a progress bar on a text terminal showing file being transferred, completion
    percentage, estimated transfer time and throughput.

    ```Progress "myfile.bin" [ 20.0 B/ 20.0 B] 100% |###############| ETA:  0:00:00 147.8 B/s```
    """

    def __init__(self, filename, filesize=None, progress_msg='Progress'):
        """
        Parameters
        ----------
        filename : :py:class:`str`
            file name to display, will be truncated for display if necessary

        filesize : :py:class:`int`
            file size in bytes

        progress_msg : :py:class:`str`
            message to display before filename, e.g. "[MyPackage] downloading", defaults to ""
        """
        if len(filename) > 34:
            filename = os.path.splitext(filename)[0][-34:]
            progress_msg = ''
        self.progress = tqdm(unit='B',
                             unit_scale=True,
                             unit_divisor=1024,
                             miniters=1,
                             total=filesize,
                             postfix={progress_msg: filename})

    def __call__(self, bytes_amount):
        self.progress.update(bytes_amount)

_GCS_BUCKET = os.environ.get("GENOMEKIT_GCS_BUCKET", "genomekit-public-dg")

def _hashfile(afile, hasher, blocksize=65536):
    """Memory efficient file hashing function.

    Parameters
    ----------
    afile : :py:class:`file`
        An open file handle to read from.

    hasher : :py:class:`_hashlib.HASH`
        A hash function, e.g. hashlib.sha256().

    blocksize : :py:class:`int`
        The size of each hash update (speed only).

    Returns
    -------
    :py:class:`str`
        The hash string.
    """
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    return hasher.digest()

class DataManager(ABC):
    def __init__(self, data_dir: str):
        """
        Implementations should accept a data directory.
        Args:
            data_dir: location where local files are cached
        """
        os.makedirs(data_dir, exist_ok=True)
        self.data_dir = data_dir

    def get_file(self, filename: str) -> str:
        """Get path to a managed data file and download if necessary.

        Parameters
        ----------
        filename
            A single filename, without a path.

        Returns
        -------
            The path to the local version of the file.
        """
        raise NotImplementedError("Descendant classes must implement get_file")

    def upload_file(self, filepath: str, filename: str, metadata: Dict[str, str]=None):
        """Upload a file to make it available in GenomeKit.

        Parameters
        ----------
        filepath
            Full path to the file to be uploaded.

        filename
            The remote filename.

        metadata
            Extra metadata to attach to the file.
        """
        raise NotImplementedError("Descendant classes must implement upload_file")

    def list_available_genomes(self):
        """List all available genomes in the data manager"""
        raise NotImplementedError("Descendant classes must implement list_available_genomes")


def _remote_equal(blob: storage.Blob, file_path: Path) -> bool:
    remote_checksum = base64.b64decode(blob.md5_hash)
    with open(file_path, 'rb') as f:
        local_checksum = _hashfile(f, hashlib.md5())
    return remote_checksum == local_checksum


class CallbackIOWrapper(ObjectWrapper, object):
    def __init__(self, callback, stream):
        super(CallbackIOWrapper, self).__init__(stream)

        write = getattr(stream, "write")

        @wraps(write)
        def wrapped_write(data, *args, **kwargs):
            res = write(data, *args, **kwargs)
            callback(data)
            return res

        self.wrapper_setattr("write", wrapped_write)

        read = getattr(stream, "read")

        @wraps(read)
        def wrapped_read(*args, **kwargs):
            data = read(*args, **kwargs)
            callback(data)
            return data

        self.wrapper_setattr("read", wrapped_read)

@contextmanager
def FileIO(filename, desc, mode, size, quiet):
    with open(filename, mode) as fp:
        with tqdm(
                disable=quiet,
                total=size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                miniters=1,
                desc=desc,
        ) as bar:
            decorated = CallbackIOWrapper(
                lambda x: bar.update(len(x)), fp
            )
            yield decorated

class DefaultDataManager(DataManager):
    """A minimal data manager implementation that retrieves files from the GCS samples bucket.
    Uploads are not supported"""
    def __init__(self, data_dir: str):
        """
        Args:
            data_dir: location where local files are cached
        """
        super().__init__(data_dir)

    @property
    def bucket(self):
        if not hasattr(self, "_bucket"):
            gcloud_client = storage.Client()
            try:
                self._bucket = gcloud_client.bucket(_GCS_BUCKET, user_project=os.environ.get("GENOMEKIT_GCS_BILLING_PROJECT", None))
            except Exception as e:
                # give the user a hint in case of permission errors
                print(e, file=sys.stderr)
                raise

        return self._bucket

    def get_file(self, filename: str) -> str:
        local_path = Path(self.data_dir, filename)

        if local_path.exists():
            return str(local_path)

        try:
            blob = self.bucket.blob(filename)
            if not blob.exists():
                raise FileNotFoundError(f"File '{filename}' not found in the GCS bucket")
        except Exception as e:
            # give the user a hint in case of permission errors
            print(e)
            raise

        # form a temporary filename to make the download safe
        temp_file = tempfile.NamedTemporaryFile(delete=False, mode="wb", dir=self.data_dir, prefix=filename, suffix=".part")
        temp_file_path = str(Path(self.data_dir, temp_file.name))
        try:
            temp_file.close()
            with FileIO(temp_file_path, filename, "wb", blob.size, quiet=False) as f:
                blob.download_to_file(f)
        except:
            os.remove(temp_file_path)
            raise
        # atomically (on POSIX) rename the file to the real one
        os.rename(temp_file_path, local_path)

        return str(local_path)

    def upload_file(self, filepath: str, filename: str, metadata: Dict[str, str]=None):
        blob = self.bucket.blob(filename)

        if blob.exists():
            blob.reload()
            if _remote_equal(blob, Path(filepath)):
                logger.info(f"File '{filename}' already exists in the GCS bucket and is identical, skipping.")
                return

            if "GK_UPLOAD_ALLOW_OVERWRITE" not in os.environ:
                raise ValueError(f"File '{filename}' already exists in the GCS bucket."
                                 "Set GK_UPLOAD_ALLOW_OVERWRITE=1 to overwrite.")
            else:
                logger.warning(f"Overwriting '{filename}' on the server.")

            blob.metadata.update(metadata or {})
        else:
            blob.metadata = metadata

        with FileIO(filepath, filename, "rb", os.path.getsize(filepath), quiet=False) as f:
            blob.upload_from_file(f)

    @lru_cache
    def list_available_genomes(self):
        blobs = self.bucket.list_blobs(match_glob="*.{2bit,cfg}")
        return [blob.name.rpartition(".")[0] for blob in blobs]
