# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import os

##################################################################


def makedirs(path):  # pragma: no cover
    """
    Calls os.makedirs, but only if the directory doesn't
    already exist. Returns the path string.

    Since makedirs can intermittently fail if the path
    directory is being scanned by Dropbox or BTSync,
    this function retries multiple times.
    """
    trials = 0
    while not os.path.exists(path):
        try:
            os.makedirs(path)
            break
        except OSError:
            trials += 1
            if trials == 4:
                raise
    return path


##################################################################


# Copied from BioPython
def _maketrans(complement_mapping):
    """Makes a python string translation table (PRIVATE).
    Arguments:
        - complement_mapping - a dictionary such as ambiguous_dna_complement
          and ambiguous_rna_complement from Data.IUPACData.
    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.
    Compatible with lower case and upper case sequences.
    For internal use only.
    """
    before = ''.join(complement_mapping.keys())
    after = ''.join(complement_mapping.values())
    before += before.lower()
    after += after.lower()
    return str.maketrans(before, after)


_DNA_COMPLEMENT_TABLE = _maketrans({
    'A': 'T',
    'B': 'V',
    'C': 'G',
    'D': 'H',
    'G': 'C',
    'H': 'D',
    'K': 'M',
    'M': 'K',
    'N': 'N',
    'R': 'Y',
    'S': 'S',
    'T': 'A',
    'V': 'B',
    'W': 'W',
    'X': 'X',
    'Y': 'R'
})


def reverse_complement(seq):
    """Computes the reverse complement of a sequence (internal)."""
    return str(seq[::-1]).translate(_DNA_COMPLEMENT_TABLE)
