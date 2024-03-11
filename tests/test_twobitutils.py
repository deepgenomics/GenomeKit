# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import platform
import unittest
import os
from tempfile import mkstemp

import pytest

import genome_kit._twobitutils as _twobitutils
if platform.system() != 'Windows':
    from twobitreader import TwoBitFile


def _check_write(test, seqs):
    # Create named temp file;  immediately close handle so can be re-opened
    fd, tmpfile = mkstemp()
    os.close(fd)
    try:
        # Dump the sequences to a 2bit file
        _twobitutils.write2bit(tmpfile, seqs)
    finally:
        os.unlink(tmpfile)  # Remove temporary file

@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
def _check_write_read(test, seqs):
    # Create named temp file;  immediately close handle so can be re-opened
    fd, tmpfile = mkstemp()
    os.close(fd)
    twobit = None
    try:
        # Dump the sequences to a 2bit file
        _twobitutils.write2bit(tmpfile, seqs)

        # Read the sequences back and ensure twobitreader extracts the original sequences exactly
        twobit = TwoBitFile(tmpfile)
        test.assertEqual(set(seqs.keys()), set(twobit))  # Same set of sequence names
        for name, seq in seqs.items():
            tbseq = twobit[name]
            test.assertEqual(len(tbseq), len(seq))
            test.assertEqual(tbseq[:].upper(), seq.upper())
            test.assertEqual(tbseq[:], seq)

    finally:
        if twobit is not None:
            # Get TwoBitFile to release the handle to temporary file
            twobit._file_handle.close()  # ("del twobit" doesn't work)
            os.unlink(tmpfile)  # Remove temporary file


@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
class TwoBitWriteEmptyTest(unittest.TestCase):
    def test(self):
        seqs = {'x': ''}
        _check_write_read(self, seqs)

        seqs = {'x': '', 'y': ''}
        _check_write_read(self, seqs)


@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
class TwoBitWriteSimpleTest(unittest.TestCase):
    def test(self):
        # Multiple of 4
        seqs = {'chr1': 'ACGTAACCGGTT', 'chr2': ''}
        _check_write_read(self, seqs)

        # Not multiple of 4
        seqs = {'chr1': 'ACGTAACCGGTTA'}
        _check_write_read(self, seqs)


@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
class TwoBitWriteMaskTest(unittest.TestCase):
    def test(self):
        seqs = {'chr1': 'ACgtAACCggTTA'}
        _check_write_read(self, seqs)

        seqs = {'chr1': 'acGTaaccGGtta'}
        _check_write_read(self, seqs)


@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
class TwoBitWriteNBlockTest(unittest.TestCase):
    def test(self):
        seqs = {'chr1': 'NNGTNNNNGGNNN'}
        _check_write_read(self, seqs)

        seqs = {'chr1': 'ACNNAACCNNTTA'}
        _check_write_read(self, seqs)


@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
class TwoBitWriteMaskNBlockTest(unittest.TestCase):
    def test(self):
        seqs = {'chr1': 'ACannNACNnngTT'}
        _check_write_read(self, seqs)

        seqs = {'chr1': 'NNncgANNTccnNN'}
        _check_write_read(self, seqs)


@pytest.mark.skipif(platform.system() == 'Windows', reason="twobitreader not available on Windows")
class TwoBitWriteMultiSequenceTest(unittest.TestCase):
    def test(self):
        seqs = {'chr1': 'ACannNACNnngTT', 'chr1_extra': 'TTgcCNCGGGGNNnnnnccaAA'}
        _check_write_read(self, seqs)


class TwoBitWriteInvalidSequenceTest(unittest.TestCase):
    def test(self):
        seqs = {'chr1': 'AAAAAAPPPAAAAAAAAA'}
        with self.assertRaises(Exception):
            _check_write(self, seqs)

        seqs = {'chr1': 'aaaaaapppaaaaaaaaa'}
        with self.assertRaises(Exception):
            _check_write(self, seqs)


if __name__ == "__main__":
    unittest.main()
