# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
#
# This file allows the genome_kit package to be run as a script.
# The script provides subcommands related to GenomeKit development
# and maintenance, such as building internal data files.
#
from __future__ import annotations

import os
import pickle
import re
import tarfile
from collections import defaultdict
from contextlib import suppress
from tempfile import TemporaryDirectory

from . import Genome, _build_annotations, gk_data
from ._gk_data_config import get_appris_filename

# Original sources from
# http://apprisws.bioinfo.cnio.es/pub/releases/2024_10.v49
# http://apprisws.bioinfo.cnio.es/pub/releases/2023_08.v48
# http://apprisws.bioinfo.cnio.es/pub/releases/2022_07.v47
# http://apprisws.bioinfo.cnio.es/pub/releases/2019_09.v29
# http://apprisws.bioinfo.cnio.es/pub/releases/2018_12.v28
# http://apprisws.bioinfo.cnio.es/pub/releases/2018_02.v27
# http://apprisws.bioinfo.cnio.es/pub/releases/2017_06.v23
# Steps:
#     - download the remote file, read it into memory
#     - create a dictionary to index the APPRIS data by geneID
#     - create a new column with the name `principalityIndex` which contains integers that range 0..6 based
#       on the principality string
#     - create a dict to map transcript ids to indices in GK's transcript table, and
#       create a dict to map gene ids to indices in GK's gene table
#     - create [[appris_indices], [appris_transcripts], [appris_transcripts_by_gene]]
#     - save the list as a .pkl file
#     - upload .pkl file to with the DataManager

# Use path of _build_appris.py
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'tests', 'data')
assert os.path.isdir(TEST_DATA_DIR)


def get_appris_principal_path(annotation):
    # see http://appris.bioinfo.cnio.es/#/help/species

    MISSING = (None, None, None)

    def parse_ensembl(name):
        if name == "Rnor_6.0.88":
            return ("rattus_norvegicus", "e88v22", "2017_06.v23")
        elif name == "Sscrofa11.1.98":
            return ("sus_scrofa", "e98v29", "2019_09.v29")
        return MISSING

    def parse_gencode(name):
        name = name[1:]
        species = "mus_musculus" if name.startswith("M") else "homo_sapiens"

        if name.startswith("19"):
            return (species, "g19v22", "2017_06.v23")
        elif name.startswith("25"):
            return (species, "e87v22", "2017_06.v23")
        elif any(name.startswith(x) for x in ["M15", "26"]):
            return (species, "e88v22", "2017_06.v23")
        elif name.startswith("27"):
            return (species, "e91v27", "2018_02.v27")
        elif any(name.startswith(x) for x in ["M19", "29"]):
            return (species, "e94v28", "2018_12.v28")
        elif any(name.startswith(x) for x in ["M30", "41"]):
            return (species, "e107v47", "2022_07.v47")
        elif name.startswith("46"):
            return (species, "e112v48", "2024_10.v49")
        elif any(name.startswith(x) for x in ["M36", "47"]):
            return (species, "e113v49", "2024_10.v49")
        return MISSING

    def parse_ncbi(name):
        if name == "v105.20190906":
            return ("homo_sapiens", "rs105v22", "2017_06.v23")
        elif name == "v109":
            return ("homo_sapiens", "rs109v28", "2018_12.v28")
        elif name == "m38.v106":
            return ("mus_musculus", "rs106v29", "2019_09.v29")
        elif name == "v110":
            return ("homo_sapiens", "rs110v48", "2023_08.v48")
        return MISSING

    source, _, remainder = annotation.partition(".")

    species, data_version, appris_version = {
        "ensembl": parse_ensembl,
        "gencode": parse_gencode,
        "ncbi_refseq": parse_ncbi,
        "ucsc_refseq": lambda _: ("homo_sapiens", "rs105v22", "2017_06.v23"),
    }.get(source, lambda _: MISSING)(remainder)
    if species is None:
        raise ValueError(f"No APPRIS defined for: {annotation}")
    # TODO: handle appris version > 99
    assert int(appris_version[-2:]) >= int(data_version[-2:])

    return (
        "apprisws.bioinfo.cnio.es/pub/releases"
        f"/{appris_version}/datafiles/{species}/{data_version}"
        "/appris_data.principal.txt")


APPRIS_DATA_HEADER = ['gene', 'geneID', 'transcriptID', 'ccds', 'principality']
OUTPUT_DATA_HEADER = ['geneID', 'transcriptID', 'principalityIndex']

# archived via: wget -q -r -np -A appris_data.principal.txt --accept-regex 'homo_sapiens|mus_musculus|rattus_norvegicus|sus_scrofa' https://apprisws.bioinfo.cnio.es/pub/releases/2024_10.v49/datafiles/
# and merged with previous backups.
APPRIS_ARCHIVE_NAME = "apprisws.bioinfo.cnio.es.2024_10.v49.tar.gz"

ISO_DICT = {
    'PRINCIPAL:1': 0,
    'PRINCIPAL:2': 1,
    'PRINCIPAL:3': 2,
    'PRINCIPAL:4': 3,
    'PRINCIPAL:5': 4,
    'ALTERNATIVE:1': 5,
    'ALTERNATIVE:2': 6
}


def _download_appris_data(url, delimiter='\t'):
    """Fetches the APPRIS data from a server (if needed), breaks down each line
    into a list, and returns the line by line APPRIS data in the form of a list.

    Parameters
    ----------
    url : str
        The URL to a `appris_data.principal.txt` file. The APPRIS data is a 5 column tab
        separated file, where columns are [geneName, geneID, transcriptID, CDS, principality]
    delimiter : str
        Optional. The tabular separator for the source APPRIS data
    """
    # Pull the file rom original source, if needed, and parse the file.
    with open(url, "r") as f:
        data = [line.rstrip().split(delimiter) for line in f]
        if data[0][0] == "Gene name":
            data = data[1:]
    return data


def _unversioned_id(id: str) -> str:
    return id.partition('.')[0]


def _process_appris_data(data):
    """Parses the principalityNum from the data.

    Parameters
    ----------
    data : list
        A list of lists where each list is in the form of [geneName, geneID, transcriptID, CDS, principality]

    Returns
    -------
    dict
        A dictionary of dictionaries in the form of
            {
                "geneID": {"transcriptID": principalityNum, ...},
                ...
            }
    """

    # Filter data and add a processed transcript_id for keying result:
    # [gene_id, transcript_id, transcript_id_num, principality_num]
    # skip MANE lines (PRINCIPAL:M, ALTERNATIVE:M)
    data = ((line[1], line[2], ISO_DICT[line[4]]) for line in data if line[4] in ISO_DICT)

    # Take data in the form [[gene_id, transcript_id, transcript_id_num, principality_num] ... ]
    result = defaultdict(dict)
    for gene_id, transcript_id, principality_num in data:
        transcript_id = _unversioned_id(transcript_id)
        gene_id = _unversioned_id(gene_id)
        # ensure there are no duplicate transcript_ids per gene
        old = result[gene_id].setdefault(transcript_id, principality_num)
        assert old == principality_num, f"duplicate transcript_id {transcript_id}"

    return result


def _get_appris_data(url, delimiter='\t'):
    """Shorthand for both downloading and processing appris data URL."""
    return _process_appris_data(_download_appris_data(url, delimiter=delimiter))


def _gen_gk_appris_data(data, genome):
    """Generates the data dependency that GK needs for fast retrieval of APPRIS scores.

    Parameters
    ----------
    data : dict
        A dictionary of geneIDs mapped to principalities in the form:
        {
            "geneID": {"transcriptID": principality, ...},
            ...
        }
    genome : :py:class:~genome_kit:Genome
        A genome object with the matching annotations for the APPRIS data

    Returns
    -------
    tuple
        [0] a list of APPRIS principality scores where each score is mapped to a
            :py:class:~genome_kit.Transcript in :py:class:~genome_kit.TranscriptTable
        [1] a list of indices for :py:class:~genome_kit.Transcript objects for fast access to all
            the transcripts with available APPRIS scores
        [2] a list of lists, where each list (ordered by APPRIS principality scores)
            of :py:class:~genome_kit.Transcript objects indices is mapped to a
            :py:class:~genome_kit.Gene object index.
    """

    # https://www.ensembl.org/info/genome/stable_ids/
    ensembl_transcript_prefix_pattern = r"ENS(\w\w\w)?TR?"

    def refseq(transcriptID):
        return _unversioned_id(transcriptID).split('_')[1]

    def gencode(transcriptID):
        return re.sub(ensembl_transcript_prefix_pattern, "", transcriptID)

    # Determine if the APPRIS data is for RefSeq or Gencode based on the first line
    transcript_id = next(iter(next(iter(data.values()))))
    if re.fullmatch(ensembl_transcript_prefix_pattern + r'\d+', transcript_id):
        simplify_transcript = gencode
    elif re.fullmatch(r'[A-Z]{2}_\d+', transcript_id):
        simplify_transcript = refseq

    gene_dict = defaultdict(list)

    # Map geneID to its index in GeneTable
    for i, gene in enumerate(genome.genes):
        gid = _unversioned_id(gene.id)
        gene_dict[gid].append(i)

    # Desired products used by genomekit
    # A list of APPRIS principality scores
    appris_indices = [None] * len(genome.transcripts)
    appris_transcripts = []
    appris_transcripts_by_gene = [[] for i in range(len(genome.genes))]

    for gene, gene_indices in gene_dict.items():
        try:
            _appris_transcripts = data[gene]
        except KeyError:
            # gene in GK is not found in APPRIS
            continue

        for gidx in gene_indices:
            for t in genome.genes[gidx].transcripts:
                try:
                    _p = _appris_transcripts[_unversioned_id(t.id)]
                except KeyError:
                    continue
                tidx = genome.transcripts.index_of(t)
                appris_indices[tidx] = _p
                appris_transcripts.append(tidx)
                appris_transcripts_by_gene[gidx].append(tidx)
            appris_transcripts_by_gene[gidx].sort(
                key=lambda tidx: (
                    appris_indices[tidx],
                    simplify_transcript(genome.transcripts[tidx].id),
                )
            )

    return appris_indices, appris_transcripts, appris_transcripts_by_gene


def build_full_appris_file(annotation, source_file, upload: bool = False):
    """Builds the full GenomeKit APPRIS binary data dependencies for a particular
    annotation and optionally uploads to the store.

    The file will be placed in the ``gk_data._config['DATA_DIR']`` directory.

    Parameters
    ----------
    annotation : str
        A GenomeKit supported genome annotation type, e.g. "gencode.v29"
    upload:
        Optional. Whether to upload to the store or not
    """

    # Path to final output file
    output_filename = get_appris_filename(annotation)
    output_filepath = os.path.join(gk_data._config['DATA_DIR'], output_filename)

    # Get APPRIS data as dict structure
    try:
        data = _get_appris_data(source_file)
    except:  # noqa
        print('Error building {}'.format(annotation))
        raise

    # Build APPRIS structures specific to this GenomeKit annotation
    genome = Genome(annotation)
    gkresults = _gen_gk_appris_data(data, genome)

    # Dump the final annotation-specific structures to disk.
    with open(output_filepath, 'wb') as fp:
        pickle.dump(gkresults, fp, 2)

    print("Wrote", output_filepath, end="")

    if upload is True:
        print(" uploading ...", end="")
        gk_data.upload_file(output_filepath, output_filename)

    print()


def build_full_appris_files(upload: bool = False):
    archive = gk_data.get_file(APPRIS_ARCHIVE_NAME)
    tar = tarfile.open(archive)

    with TemporaryDirectory() as tmpdir:
        tar.extractall(tmpdir)

        for annotation in gk_data.data_manager.list_available_genomes():
            with suppress(ValueError):
                build_full_appris_file(annotation, os.path.join(tmpdir, get_appris_principal_path(annotation)), upload)


def build_test_appris_file(annotation, source_file):
    """Builds test GenomeKit APPRIS binary files, using mini-annotations for testing."""

    # Create a Genome object to load mini-annotations
    data_dir = os.path.join(TEST_DATA_DIR, "mini1")

    data = _get_appris_data(source_file)
    genome = Genome(annotation)
    output = _gen_gk_appris_data(data, genome)
    output_filepath = os.path.join(data_dir, get_appris_filename(annotation))

    with open(output_filepath, 'wb') as fp:
        pickle.dump(output, fp, 2)

    print("Wrote", output_filepath)


def build_test_appris_files():
    all_test_annotations = list(
        _build_annotations.GENCODE_OR_NCBI_TEST_ANNOTATIONS.values()
    ) + list(_build_annotations.UCSC_REFSEQ_TEST_ANNOTATIONS.values())

    archive = gk_data.get_file(APPRIS_ARCHIVE_NAME)
    tar = tarfile.open(archive)

    with TemporaryDirectory() as tmpdir:
        tar.extractall(tmpdir)

        for annotation in all_test_annotations:
            if annotation in ["gencode.v29lift37", "ncbi_refseq.hg38.p14_RS_2024_08"]:
                continue  # leave this out to test missing APPRIS (test_appris_transcripts_unavailable)
            build_test_appris_file(annotation + ".mini" , os.path.join(tmpdir, get_appris_principal_path(annotation)))
