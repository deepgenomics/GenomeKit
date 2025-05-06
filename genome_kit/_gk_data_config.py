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


def get_appris_filename(annotation):
    return f"appris.{annotation}.v{APPRIS_BINARY_VERSION}.pkl"


# MANE data file names
# 1 - preprocessed MANE data
#
# NOTE: if you bump this, then `python -m genome_kit build --mane --test-mane`
#       and `git rm` the old test files
MANE_BINARY_VERSION = 1


def get_mane_filename(annotation):
    return (
        f"mane.{annotation}.v{MANE_BINARY_VERSION}.pkl"
    )
