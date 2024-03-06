# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import os
from pathlib import Path

# TODO: The C++ backend code currently uses GENOMEKIT_QUIET environment variable
#       to suppress printing progress messages. Would be good to replace this
#       with proper logging integration.
os.environ["GENOMEKIT_QUIET"] = "1"

from genome_kit import Genome, GenomeDNA, gk_data

test_data_dir = Path(__file__).parent / "data" / "mini1"

_get_file_original = gk_data.data_manager.get_file


def get_test_file(filename):
    if (test_data_dir / filename).exists():
        return str(test_data_dir / filename)
    return _get_file_original(filename)


gk_data.data_manager.get_file = get_test_file


#####################################################


class MiniGenomeDNA(GenomeDNA):
    """A GenomeDNA object accessing a small test file, specified by name."""

    def __init__(self, name=None):
        return super(MiniGenomeDNA, self).__init__(str(test_data_dir / f"{name}.2bit"))


#####################################################


class MiniGenome(Genome):
    """A Genome object accessing small test files, specified by name."""

    def __init__(self, config=""):
        return Genome.__init__(self, config)

    def __new__(cls, config=""):
        if not config:
            config = 'gencode.v29lift37'  # hg19 + gencode by default

        config = f"{config}.mini" if 'test_genome' not in config and not config.endswith(".mini") else config

        return Genome.__new__(cls, config)


#####################################################


def check_pythonic_indexing(test, obj):
    """Tests that `obj` supports Python-style sequence indexing.
    This method would normally be called from a unit test method,
    so pass `self` as `test`."""
    n = len(obj)
    assert n > 0, "Can't test empty sequence"
    indices = sorted(set(list(range(min(n, 5))) + list(range(max(0, n - 5), n))))
    for i in indices:
        test.assertEqual(obj[-i - 1], obj[n - i - 1])
    with test.assertRaises(IndexError):
        obj[n]
    with test.assertRaises(IndexError):
        obj[-n - 1]


#####################################################


def dumptext(outfile, *texts):
    """Dumps a each text string to a file.

    Leading newline/whitespace is stripped from each line, so that tests can
    indent their Python string literals to match surround code, whilst
    having that indentation stripped from the resulting file.
    """
    with open(outfile, "w") as f:
        for text in texts:
            text = text.lstrip()
            text = "\n".join([line.lstrip(" ") for line in text.split("\n")])
            f.write(text)
    return outfile
