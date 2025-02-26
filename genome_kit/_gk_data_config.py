# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import os

import appdirs

_config = {
    # DATA_DIR is taken either from the GENOMEKIT_DATA_DIR environment variable
    # or determined in a platform-dependent way using the appdirs package.
    # E.g. on a Mac, DATA_DIR will be ~/Library/Application Support/genome_kit.
    # Multiple directories can be specified via a colon seperated list.
    "DATA_DIR": os.environ.get(
        "GENOMEKIT_DATA_DIR", appdirs.user_data_dir("genome_kit")
    ),
}

#: Utility function for adding several managed file names all based on the same format.
# APPRIS data file names
# 1 - preprocessed APPRIS data
# 2 - build from APPRIS source
# 3 - NCBI refseq v109 uses type==pseudogene for some genes, including them caused indices to shift
# 4 - NCBI use refseq sources (which are archived) instead of all(GCF), but some indices are swapped
# 5 - Differentiate genes/transcripts in UCSC mapping to multiple regions
# 6 - Fix accidental removal of sorting by principality
# 7 - Allow all contigs
#
# NOTE: if you bump this, then `python -m genome_kit build --appris --test-appris`
#       and `git rm` the old test files
APPRIS_BINARY_VERSION = 7


def get_appris_version(annotation):
    if any(x in annotation for x in ["Sscrofa11.1.98", ".m38.v"]):
        return "2019_09.v29"
    elif any(
        x in annotation for x in ["gencode.v29", "gencode.vM19", "ncbi_refseq.v109"]
    ):
        return "2018_12.v28"
    elif any(x in annotation for x in ["gencode.v41", "gencode.vM30"]):
        return "2022_07.v47"
    elif any(x in annotation for x in ["ncbi_refseq.v110"]):
        return "2023_08.v48"
    elif any(x in annotation for x in ["gencode.v27"]):
        return "2018_02.v27"
    return "2017_06.v23"


def get_appris_filename(annotation):
    return f"appris.{get_appris_version(annotation)}_{annotation}.v{APPRIS_BINARY_VERSION}.pkl"


# MANE data file names
# 1 - preprocessed MANE data
#
# NOTE: if you bump this, then `python -m genome_kit build --mane --test-mane`
#       and `git rm` the old test files
MANE_BINARY_VERSION = 1

_SUPPORTED_MANE_VERSIONS_BY_ANNO = {
    "ncbi_refseq.hg38.p14_RS_2024_08": "1.4",
    "ncbi_refseq.hg38.p14_RS_2024_08.mini": "1.4",
    "gencode.v41": "1.0",
    "gencode.v41.mini": "1.0",
}


def get_mane_version(annotation):
    """
    Get the MANE version based on the annotation version.
    :raises ValueError: if there is no matching MANE version for the specified annotation
    """
    try:
        return _SUPPORTED_MANE_VERSIONS_BY_ANNO[annotation]
    except KeyError:
        raise ValueError(f"MANE not supported for annotation {annotation}")


def get_mane_filename(annotation):
    return (
        f"mane.{get_mane_version(annotation)}_{annotation}.v{MANE_BINARY_VERSION}.pkl"
    )
