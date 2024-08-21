# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import gc
from tempfile import mkstemp
from collections import Counter
from itertools import product

import numpy as np
import unittest

from genome_kit import GenomeTrack
from genome_kit import GenomeTrackBuilder
from genome_kit import Interval

from . import MiniGenome
from . import dumptext


# Generate a data block of numbers for a given interval
def make_data(interval, dim, pattern, dtype):
    # Make an array of the right size
    s, e = interval.start, interval.end
    r = np.empty((e - s, dim), dtype=dtype)

    # Fill each column with numbers that depend on both their
    # position and on which track dimension they correspond to.
    for d in range(dim):
        r[:, d] = pattern[np.arange(s + d, e + d) % len(pattern)]
    return r


def as_list(array):
    return array.reshape(array.shape[0]).tolist()


def make_builder(outfile, etype="f16", strandedness="single_stranded", reference_genome=MiniGenome("hg19"), **kwargs):
    return GenomeTrackBuilder(outfile, etype, strandedness, reference_genome, **kwargs)

########################################################


class TestBuildTrack(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        # Create a tmpfile that the test will continually clobber as it builds various tracks
        fd, cls.tmpfile = mkstemp()
        os.close(fd)
        fd, cls.tmpfile2 = mkstemp()
        os.close(fd)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.tmpfile)  # Remove temporary file
        os.unlink(cls.tmpfile2)

    def make_builder(self, etype="f16", strandedness="single_stranded", reference_genome=MiniGenome("hg19"), **kwargs):
        return make_builder(self.tmpfile, etype=etype, strandedness=strandedness, reference_genome=reference_genome, **kwargs)

    def check_track(self, interval, expected_values, dtype=None):
        with GenomeTrack(self.tmpfile) as t:
            if dtype:
                values = t(interval, dtype)
            else:
                values = t(interval)

            np.testing.assert_equal(expected_values, values)

    def test_errors(self):
        # GenomeTrackBuilder
        with self.assertRaises(TypeError):
            make_builder(None)  # Invalid file path
        with self.assertRaises(ValueError):
            self.make_builder(dim=0)  # Invalid dim
        with self.assertRaises(ValueError):
            self.make_builder(dim=-1)  # Invalid dim
        with self.assertRaises(TypeError):
            self.make_builder(None)  # Invalid etype
        with self.assertRaises(ValueError):
            self.make_builder('foo')  # Unrecognized etype
        with self.assertRaises(ValueError):
            self.make_builder('m0', dim=2)  # Can't create 2D mask track
        with self.assertRaises(TypeError):
            self.make_builder(reference_genome='hg19')  # str instead of Genome
        track = self.make_builder(dim=2)
        with self.assertRaisesRegex(ValueError, "stride"):
            track.set_data(Interval('chr1', '+', 0, 20, track.refg), np.zeros((20, 2), np.float16)[::-1])
        with self.assertRaisesRegex(ValueError, "stride"):
            track.set_data(Interval('chr1', '+', 0, 20, track.refg), np.zeros((20, 2), np.float16)[:, ::-1])

        # Build a valid .gtrack file for some tests
        builder = self.make_builder()
        with self.assertRaises(ValueError):
            builder.set_data(Interval('chr1', '+', 0, 20, builder.refg), np.zeros((10, 1), np.float16))
        with self.assertRaises(ValueError):
            builder.set_data(Interval('chr1', '+', 0, 20, builder.refg), np.zeros((30, 1), np.float16))
        builder.set_data(Interval('chr1', '+', 0, 20, builder.refg), np.zeros((20, 1), np.float16))
        builder.finalize()

        # GenomeTrack
        with self.assertRaises(ValueError):
            GenomeTrack(None)  # Invalid file path

        # Should be able to instantiate track with invalid file name,
        # but raise error when trying to access the track (open on demand)
        track = GenomeTrack("file_that_does_not_exist.gtrack")  # Invalid file path
        with self.assertRaises(IOError):
            track(Interval('chr1', '+', 0, 20, track.refg))

        # TODO: check mask special cases, ensure they throw errors
        #

        # TODO: ensure encoding value out of range throws error
        #       for float8, float8u, float16

        # TODO: ensure set_data_from_wig works for both single and double stranded
        #       (generate a temp wig file)

    def test_builder(self):
        genome = MiniGenome("hg19")

        # Intervals on both strands. The general pattern is shown
        # below, where the number represents the position modulo 10
        #
        #   [+]   0123456789-----5678901234-----
        #   [-]   -----5678901234-----0123456789
        #
        #
        all_intervals = [
            Interval('chr1', '+', 0, 10, genome),
            Interval('chr1', '+', 15, 25, genome),
            Interval('chr1', '-', 5, 15, genome),
            Interval('chr1', '-', 20, 30, genome),
        ]

        # (etype, dtype, default_value, pattern)
        # where pattern is used to draw values from when filling the track
        test_cases = [
            ('m0',   np.bool_,    False, np.array([ 1,  1,  1,  1])),
            ('u1',   np.uint8,        0, np.array([ 0,  1,  0,  1])),
            ('u2',   np.uint8,        3, np.array([ 0,  1,  2,  3])),
            ('u3',   np.uint8,        7, np.array([ 3,  4,  5,  6])),
            ('u4',   np.uint8,       15, np.array([ 5,  6,  7,  8])),
            ('u5',   np.uint8,       31, np.array([20, 21, 22, 23])),
            ('u6',   np.uint8,       63, np.array([50, 51, 52, 53])),
            ('u8',   np.uint8,      127, np.array([ 5,  6,  7,  8])),
            ('u8',   np.float32,    127, np.array([ 5,  6,  7,  8])),
            ('i8',   np.int8,      -128, np.array([-5, -3,  3,  5])),
            ('i8',   np.float32,   -128, np.array([-5, -3,  3,  5])),
            ('f2',   np.float16, -65504, np.array([1024, 1025, 1026, 1027])),  # Match fdict
            ('f3',   np.float16, -65504, np.array([1024, 1025, 1026, 1027])),
            ('f4',   np.float16, -65504, np.array([1024, 1025, 1026, 1027])),
            ('f5',   np.float16, -65504, np.array([1024, 1025, 1026, 1027])),
            ('f6',   np.float16, -65504, np.array([1024, 1025, 1026, 1027])),
            ('f8',   np.float16, -65504, np.array([1024, 1025, 1026, 1027])),
            ('f16',  np.float16, -65504, np.array([-16.5, -1.25,  1.25, 16.5])),
            ('f16',  np.float32, -65504, np.array([-16.5, -1.25,  1.25, 16.5])),
            ('f32',  np.float32,   6e6, np.array([-1e6, -2e-6, 3e-6, 4e6])),
        ]  # yapf: disable

        fdict = np.arange(1024, 1024 + 256, dtype=np.float32)

        for strandedness in ("single_stranded", "strand_unaware"):
            for etype, dtype, default_value, pattern in test_cases:
                dims = (1, 2) if etype != 'm0' else (1, )  # 'mask' only supports 1D tracks
                for dim in dims:
                    for res in (1, 2):

                        ######################################
                        # If stranded, use + and - intervals, otherwise only + intervals
                        intervals = all_intervals if strandedness != "single_stranded" else [
                            _ for _ in all_intervals if _.is_positive_strand()
                        ]

                        ######################################
                        # Create a numpy array representing the track in its entirety.
                        track_len = max([interval.end for interval in all_intervals])
                        enc_data = {'+': np.empty((track_len, dim), dtype), '-': np.empty((track_len, dim), dtype)}
                        enc_data['+'][:] = default_value
                        enc_data['-'][:] = default_value
                        for i in intervals:
                            enc_data[i.strand][i.start:i.end, :] = make_data(i, dim, pattern, dtype)
                        if strandedness == "single_stranded":
                            enc_data['-'] = enc_data['+']

                        ########################################
                        # Build a single-stranded track with 2 values per position.
                        # At the same time, fill our 'ground truth' track with values.
                        build = self.make_builder(etype, strandedness, genome, dim=dim, resolution=res)

                        # Configure the builder based on current test settings
                        build.set_default_value(default_value)

                        # Set the dictionary, if relevant
                        if etype in ('f2', 'f3', 'f4', 'f5', 'f6', 'f8'):
                            build.set_dict(fdict[:2**int(etype[1])])
                        if etype.startswith("f"):
                            build.set_clamping()

                        # Fill the builder with data
                        for i in intervals:
                            data = make_data(i, dim, pattern, dtype) if etype != 'm0' else None
                            scaled = Interval.from_dna0(i.chrom, i.strand, i.start * res, i.end * res, i.refg)
                            build.set_data(scaled, data)

                        # Dump data and index to disk
                        build.finalize()

                        ######################################
                        # Load the track and check the entire range of values.
                        with GenomeTrack(self.tmpfile) as track:
                            interval = Interval('chr1', '+', 0, track_len * res, track.refg)

                            dec_data = track(interval, dtype)
                            self.assertEqual(dec_data.dtype, dtype)
                            self.assertEqual(dec_data.shape, (track_len * res, dim))
                            np.testing.assert_equal(dec_data, np.repeat(enc_data['+'], res, axis=0))

                            dec_data = track(interval.as_negative_strand(), dtype)
                            self.assertEqual(dec_data.dtype, dtype)
                            self.assertEqual(dec_data.shape, (track_len * res, dim))
                            np.testing.assert_equal(dec_data, np.repeat(enc_data['-'], res, axis=0)[::-1])

    def test_attributes(self):
        builder = self.make_builder(dim=3)
        self.assertEqual(builder.dim, 3)
        self.assertEqual(builder.etype, 'f16')
        self.assertEqual(builder.dtype, np.float16)
        self.assertEqual(builder.filename, self.tmpfile)
        self.assertIsInstance(builder.reference_genome, str)
        self.assertIsInstance(builder.refg, str)
        self.assertIsInstance(builder.__repr__(), str)
        builder.finalize()

        self.assertIsInstance(GenomeTrack.gtrack_version(), int)

        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual(track.dim, 3)
            self.assertEqual(track.etype, 'f16')
            self.assertEqual(track.dtype, np.float16)
            self.assertEqual(track.filename, self.tmpfile)
            self.assertIsInstance(track.reference_genome, str)
            self.assertIsInstance(track.refg, str)
            self.assertIsInstance(track.__repr__(), str)

    def test_query_args(self):
        # Build a valid .gtrack file for some tests
        builder = self.make_builder()
        interval = Interval('chr1', '+', 0, 20, builder.refg)
        builder.set_data(interval, np.zeros((20, 1), np.float16))
        builder.finalize()
        track = GenomeTrack(self.tmpfile)

        # Ok, but no dtype
        self.assertEqual(track(interval).dtype, np.float16)

        # Wrong interval arg types
        with self.assertRaises(TypeError):
            track()
        with self.assertRaises(TypeError):
            track(None)
        with self.assertRaises(TypeError):
            track(20)
        with self.assertRaises(TypeError):
            track('abcde')
        with self.assertRaises(TypeError):
            track(track)

        # Positional and keyword dtype
        self.assertEqual(track(interval, np.float32).dtype, np.float32)
        self.assertEqual(track(interval, dtype=np.float32).dtype, np.float32)

        # dtype = None
        self.assertEqual(track(interval).dtype, np.float16)
        self.assertEqual(track(interval, None).dtype, np.float16)
        self.assertEqual(track(interval, dtype=None).dtype, np.float16)

        # Empty keywords
        empty_dict = {}
        self.assertEqual(track(interval, **empty_dict).dtype, np.float16)

    def test_intervals(self):
        # Build a valid .gtrack file for some tests
        builder = self.make_builder(strandedness="strand_aware")
        intervals = [Interval('chr1', '+', 10, 30, builder.refg)]
        intervals.append(intervals[0].as_opposite_strand().shift(1))
        builder.set_data(intervals[0], np.zeros((20, 1), np.float16))
        builder.set_data(intervals[1], np.ones((20, 1), np.float16))
        builder.finalize()

        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual(set(intervals), set(track.intervals))

        builder = self.make_builder(strandedness="strand_aware", resolution=5)
        intervals = [Interval('chr1', '+', 10, 30, builder.refg)]
        intervals.append(intervals[0].as_opposite_strand().shift(5))
        builder.set_data(intervals[0], np.zeros((20//5, 1), np.float16))
        builder.set_data(intervals[1], np.ones((20//5, 1), np.float16))
        builder.finalize()

        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual(set(intervals), set(track.intervals))

    def test_out(self):
        # Build a valid .gtrack file for some tests
        builder = self.make_builder(strandedness="strand_aware")
        interval = Interval('chr1', '+', 10, 30, builder.refg)
        builder.set_data(interval, np.ones((20, 1), np.float16))
        builder.finalize()

        with GenomeTrack(self.tmpfile) as track:
            out = np.zeros((20, 1), np.float16)
            res = track(interval, out=out)
            np.testing.assert_equal(out, np.ones((20, 1), np.float16))
            np.testing.assert_equal(res, np.ones((20, 1), np.float16))

            out = np.zeros((20, 1), np.float32)
            res = track(interval, out=out)
            np.testing.assert_equal(out, np.ones((20, 1), np.float32))

            out = np.zeros((40, 1), np.float16)
            res = track(interval, out=out[::2])
            np.testing.assert_equal(out[::2], np.ones((20, 1), np.float16))
            np.testing.assert_equal(out[1::2], np.zeros((20, 1), np.float16))


            with self.assertRaisesRegex(ValueError, "Dimension"):
                out = np.zeros((20, 1, 1), np.float16)
                track(interval, out=out)

            with self.assertRaisesRegex(ValueError, "Row"):
                out = np.zeros((19, 1), np.float16)
                track(interval, out=out)

            with self.assertRaisesRegex(ValueError, "Column"):
                out = np.zeros((20, 2), np.float16)
                track(interval, out=out)

            with self.assertRaisesRegex(ValueError, "writable"):
                out = np.zeros((20, 1), np.float16)
                out.flags.writeable = False
                track(interval, out=out)

            with self.assertRaisesRegex(ValueError, "Negative"):
                out = np.zeros((20, 1), np.float16)
                track(interval, out=out[::-1])

    def _try_parse_wig(self,
                       dim,
                       res,
                       wigtext,
                       data_size=None,
                       index_size=None,
                       default_value=None,
                       restriction=None,
                       clamping=False,
                       dictionary=None,
                       etype='f16',
                       reference_genome=MiniGenome("hg19")):
        dumptext(self.tmpfile2, wigtext)
        builder = self.make_builder(etype, reference_genome=reference_genome, dim=dim, resolution=res)
        refg = builder.refg
        try:
            if default_value is not None:
                builder.set_default_value(default_value)
            if restriction is not None:
                builder.set_restriction(restriction)
            if clamping is True:
                builder.set_clamping()
            if dictionary is not None:
                builder.set_dict(dictionary)
            builder.set_data_from_wig(self.tmpfile2)
            builder.finalize()
            if data_size is not None:
                self.assertEqual(builder.data_size, data_size)
            if index_size is not None:
                self.assertEqual(builder.index_size, index_size)
        finally:
            # If an exception occurred during set_data_from_wig, we may
            # not reach builder.finalize() and the builder may hold the
            # file-handle until it's collected.
            # If we don't ensure the builder is collected, subsequent tests
            # may fail to open a file handle to the common tmp file.
            del builder
            while gc.collect():
                pass
        return refg

    def test_wig_parsing(self):
        def try_parse(*args, **kwargs):
            return self._try_parse_wig(*args, **kwargs)

        # Zero as a negative
        try_parse(1, 1, """
            fixedStep chrom=chr1 start=1001 step=1 span=1
            -0
            """)
        with GenomeTrack(self.tmpfile) as track:
            ref = np.array([[0.0]])
            self.assertTrue(np.array_equal(ref, track(Interval("chr1", "+", 1000, 1001, track.refg))))

        # Zero with more than one place after decimal
        try_parse(1, 1, """
            fixedStep chrom=chr1 start=1001 step=1 span=1
            0.000
            """)
        with GenomeTrack(self.tmpfile) as track:
            ref = np.array([[0.0]])
            self.assertTrue(np.array_equal(ref, track(Interval("chr1", "+", 1000, 1001, track.refg))))

        # Malformed zero
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                0abcd
                """)

        # Malformed float
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                1.23abcd
                """)

        # Empty float string
        with self.assertRaises(ValueError):
            try_parse(
                2, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                123\t
                """)

        # Valid fixedStep with dim=2 and resolution=2
        try_parse(
            2,
            2,
            """
            fixedStep chrom=chr1 start=1001 step=2 span=2
            1\t2
            3\t4
            fixedStep chrom=chr1 start=1005 step=2 span=2
            5\t6
            7\t8
            fixedStep chrom=chr1 start=1009 step=1 span=1
            9\t10
            fixedStep chrom=chr2 start=2001 step=2 span=2
            -1\t-2
            -3\t-4
            fixedStep chrom=chr2 start=2005 step=2 span=2
            -5\t-6
            -7\t-8
            fixedStep chrom=chr2 start=2009 step=1 span=1
            -9\t-10
            """,
            data_size=10,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            ref = np.repeat(np.arange(1, 11, dtype=np.float16).reshape((5, 2)), 2, axis=0)
            self.assertTrue(np.array_equal(ref, track(Interval("chr1", "+", 1000, 1010, track.refg))))
            self.assertTrue(np.array_equal(-ref, track(Interval("chr2", "+", 2000, 2010, track.refg))))

        # Valid variableStep with dim=2 and resolution=2
        try_parse(
            2,
            2,
            """
            variableStep chrom=chr1 span=2
            1001\t1\t2
            1003\t3\t4
            variableStep chrom=chr1 span=2
            1005\t5\t6
            1007\t7\t8
            variableStep chrom=chr1 span=1
            1009\t9\t10
            variableStep chrom=chr2 span=2
            2001\t-1\t-2
            2003\t-3\t-4
            variableStep chrom=chr2 span=2
            2005\t-5\t-6
            2007\t-7\t-8
            variableStep chrom=chr2 span=1
            2009\t-9\t-10
            """,
            data_size=10,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            ref = np.repeat(np.arange(1, 11, dtype=np.float16).reshape((5, 2)), 2, axis=0)
            self.assertTrue(np.array_equal(ref, track(Interval("chr1", "+", 1000, 1010, track.refg))))
            self.assertTrue(np.array_equal(-ref, track(Interval("chr2", "+", 2000, 2010, track.refg))))

        # Unknown chrom
        with self.assertRaises(ValueError):
            try_parse(
                2, 2, """
                fixedStep chrom=chrQ start=1009 step=1 span=1
                0\t0
                """
            )

        # Unexpected line value
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                0.1
                hello
                """)

        # Invalid numeric value
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                0.1
                @0.53
                """)

        # Value not encodable as f16
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                100000.0
                """)

        # f32 track
        refg = try_parse(
            1,
            1,
            """
            fixedStep chrom=chr1 start=1001 step=1 span=1
            100000.0
            """,
            etype='f32')
        self.check_track(Interval('chr1', '+', 999, 1002, refg), [[np.nan], [100000.], [np.nan]])
        with self.assertRaisesRegex(TypeError, 'Cannot decode as float16 from encoded type f32'):
            self.check_track(Interval('chr1', '+', 999, 1002, refg), [[np.nan], [100000.], [np.nan]], dtype=np.float16)

        refg = try_parse(
            1,
            1,
            """
            fixedStep chrom=chr1 start=1001 step=1 span=1
            100000.0
            """,
            default_value=200000.,
            etype='f32')
        self.check_track(Interval('chr1', '+', 999, 1002, refg), [[200000.], [100000.], [200000.]])

        # f8 track with clamping
        maxval = 40000
        d = np.linspace(1, maxval, 256).astype(np.float16)
        refg = try_parse(1, 1,
                  """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                1000.0
                fixedStep chrom=chr1 start=1002 step=1 span=1
                10000.0
                fixedStep chrom=chr1 start=1003 step=1 span=1
                40000.0
                fixedStep chrom=chr1 start=1004 step=1 span=1
                70000.0
                  """,
                  etype="f8",
                  clamping=True,
                  dictionary=d)
        fmaxval = float(maxval)
        self.check_track(
            Interval('chr1', '+', 1000, 1004, refg),
            [[d[6]], [d[64]], [fmaxval], [fmaxval]],
            # [[  942.], [10040.], [40000.], [40000.]]
            dtype=np.float16,
        )

        # Span and step do not match
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                fixedStep chrom=chr1 start=1005 step=5 span=10
                0.1
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                fixedStep chrom=chr1 start=1005 step=10 span=5
                0.1
                """)

        # Mix fixedStep and variableStep
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                0.1
                0.2
                variableStep chrom=chr1 span=1
                2001\t0.3
                2002\t0.4
                """)

        # Start not aligned to track resolution
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                fixedStep chrom=chr1 start=1005 step=10 span=10
                0.1
                0.2
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                variableStep chrom=chr1 span=10
                1001\t0.1
                1002\t0.2
                """)

        # Last span on a chromosome is irregular (doesn't match
        # track resolution) AND has multiple data
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                fixedStep chrom=chr1 start=1001 step=10 span=10
                0.1
                0.2
                fixedStep chrom=chr1 start=1021 step=5 span=5
                0.3
                0.4
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                variableStep chrom=chr1 span=10
                1001	0.1
                1011	0.2
                variableStep chrom=chr1 span=5
                1021	0.3
                1031	0.4
                """)

        # Irregular span isn't last span on a chromosome
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                fixedStep chrom=chr1 start=1001 step=10 span=10
                0.1
                0.2
                fixedStep chrom=chr1 start=1021 step=5 span=5
                0.3
                fixedStep chrom=chr1 start=1031 step=10 span=10
                0.4
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 10, """
                variableStep chrom=chr1 span=10
                1001	0.1
                1011	0.2
                variableStep chrom=chr1 span=5
                1021	0.3
                variableStep chrom=chr1 span=10
                1031	0.4
                """)

        # Negative start
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=-1 step=1 span=1
                0.1
                """)
        with self.assertRaises(ValueError):
            try_parse(1, 1, """
                variableStep chrom=chr1 span=1
                -1\t0.1
                """)

        # Well-formed line, but with invalid numeric value
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=@1001 step=1 span=1
                0.1
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                variableStep chrom=chr1 span=1
                @1001\t0.1
                """)

        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001.5 step=1 span=1
                0.1
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                variableStep chrom=chr1 span=1
                1001.5\t0.1
                """)

        # Too few columns per line
        with self.assertRaises(ValueError):
            try_parse(
                2, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                0.1
                """)
        with self.assertRaises(ValueError):
            try_parse(2, 1, """
                variableStep chrom=chr1 span=1
                1001\t0.1
                """)

        # Too many columns per line
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                fixedStep chrom=chr1 start=1001 step=1 span=1
                0.1\t0.2
                """)
        with self.assertRaises(ValueError):
            try_parse(
                1, 1, """
                variableStep chrom=chr1 span=1
                1001\t0.1\t0.2
                """)

    def test_bedgraph_parsing(self):
        def try_parse(dim, res, bedgraphtext, data_size=None, index_size=None, restriction=None, etype='f16', clamping=False, dictionary=None):
            dumptext(self.tmpfile2, bedgraphtext)
            with self.make_builder(etype, dim=dim, resolution=res) as builder:
                if restriction:
                    builder.set_restriction(restriction)
                if clamping is True:
                    builder.set_clamping()
                if dictionary is not None:
                    builder.set_dict(dictionary)
                builder.set_data_from_bedgraph(self.tmpfile2)
                builder.finalize()
                if data_size is not None:
                    self.assertEqual(builder.data_size, data_size)
                if index_size is not None:
                    self.assertEqual(builder.index_size, index_size)
                return builder.refg

        # Valid fixedStep with dim=2 and resolution=2
        try_parse(
            1,
            2,
            """
            #bedGraph section chr1:1000-1005
            chr1\t1000\t1002\t1
            chr1\t1002\t1004\t2
            chr1\t1004\t1005\t3
            #bedGraph section chr2:2000-2005
            chr2\t2000\t2002\t1
            chr2\t2002\t2004\t2
            chr2\t2004\t2005\t3
            """,
            data_size=6,
            index_size=2)  # 6 because resolution=2
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual([1, 1, 2, 2, 3], as_list(track(Interval("chr1", "+", 1000, 1005, track.refg))))
            self.assertEqual([1, 1, 2, 2, 3], as_list(track(Interval("chr2", "+", 2000, 2005, track.refg))))

        # invalid chrom
        with self.assertRaises(ValueError):
            try_parse(1, 2,
                """
                #bedGraph section chrQ:1000-1010
                chrQ\t1000\t1010\t0
                """)

        # Track dimension not 1
        with self.assertRaises(ValueError):
            try_parse(2, 1, """chr1\t1000\t1010\t1""")
        with self.assertRaises(ValueError):
            try_parse(2, 1, """chr1\t1000\t1010\t1\t2""")

        # Not enough columns
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t1000\t1010""")

        # Too many columns
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t1000\t1010\t1\t2""")

        # Malformed start/end/data
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t1000.5\t1010\t1""")
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t1000\t1010.5\t1""")
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t1000\t1010\t@1.5""")

        # Negative start
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t-1\t10\t1""")

        # End <= start
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t10\t9""")
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t10\t10""")

        # Misaligned start
        with self.assertRaises(ValueError):
            try_parse(1, 10, """chr1\t1\t10\t1""")

        # Misaligned end when not last interval in chromosome
        with self.assertRaises(ValueError):
            try_parse(1, 10, """
                chr1\t0\t9\t1
                chr1\t10\t20\t2
                """)

        # value not encodable in f16
        with self.assertRaises(ValueError):
            try_parse(1, 1, """chr1\t8\t10\t100000""", etype='f16')

        # f32 track
        refg = try_parse(1, 1, """chr1\t8\t10\t100000""", etype='f32')
        self.check_track(Interval('chr1', '+', 7, 11, refg), [[np.nan], [100000.], [100000.], [np.nan]])

        # f8 track with clamping
        maxval = 40000
        d = np.linspace(1, maxval, 256).astype(np.float16)
        refg = try_parse(1, 1,
            """
            chr1\t1000\t1001\t1000
            chr1\t1001\t1002\t10000
            chr1\t1002\t1003\t40000
            chr1\t1003\t1004\t70000
            """,
            etype="f8",
            clamping=True,
            dictionary=d)
        fmaxval = float(maxval)
        self.check_track(
            Interval('chr1', '+', 1000, 1004, refg),
            [[d[6]], [d[64]], [fmaxval], [fmaxval]],
            # [[  942.], [10040.], [40000.], [40000.]]
            dtype=np.float16,
        )

    def test_bed_parsing(self):
        def try_parse(res, strandedness, categories, bedtext, etype=None, data_size=None, index_size=None,
                      restriction=None, clamping=False, dictionary=None):
            dumptext(self.tmpfile2, bedtext)
            if etype is None:
                etype = 'f16' if categories is None else 'u4'
            with self.make_builder(etype, strandedness, resolution=res) as builder:
                builder.set_default_value(0)
                if restriction:
                    builder.set_restriction(restriction)
                if clamping:
                    builder.set_clamping()
                if dictionary is not None:
                    builder.set_dict(dictionary)
                if categories is None:
                    builder.set_data_from_bed(self.tmpfile2)
                else:
                    builder.set_data_from_bed(self.tmpfile2, categories=categories)
                builder.finalize()
                if data_size is not None:
                    self.assertEqual(builder.data_size, data_size)
                if index_size is not None:
                    self.assertEqual(builder.index_size, index_size)
                return builder.refg

        # Valid BED (resolution=2, unstranded, categories) WITHOUT strand specified
        try_parse(
            2,
            "single_stranded", ["A", "B", "C"],
            """
            chr1\t1000\t1002\tA
            chr1\t1002\t1004\tB
            chr1\t1004\t1007\tC
            chr2\t2000\t2002\tA
            chr2\t2002\t2004\tB
            chr2\t2004\t2007\tC
            """,
            data_size=8,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual([0, 0, 1, 1, 2, 2, 2], as_list(track(Interval("chr1", "+", 1000, 1007, track.refg))))
            self.assertEqual([0, 0, 1, 1, 2, 2, 2], as_list(track(Interval("chr2", "+", 2000, 2007, track.refg))))

        # Valid BED (resolution=2, unstranded, categories) WITH strand specified
        try_parse(
            2,
            "single_stranded", ["A", "B", "C"],
            """
            chr1\t1000\t1002\tA\t0\t.
            chr1\t1002\t1004\tB\t0\t.
            chr1\t1004\t1007\tC\t0\t.
            chr2\t2000\t2002\tA\t0\t.
            chr2\t2002\t2004\tB\t0\t.
            chr2\t2004\t2007\tC\t0\t.
            """,
            data_size=8,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual([0, 0, 1, 1, 2, 2, 2], as_list(track(Interval("chr1", "+", 1000, 1007, track.refg))))
            self.assertEqual([0, 0, 1, 1, 2, 2, 2], as_list(track(Interval("chr2", "+", 2000, 2007, track.refg))))

        # Valid BED (resolution=2, stranded, categories)
        try_parse(
            2,
            "strand_unaware", ["A", "B", "C"],
            """
            chr1\t1000\t1002\tA\t0\t+
            chr1\t1002\t1004\tB\t0\t+
            chr1\t1004\t1007\tC\t0\t+
            chr1\t1002\t1004\tA\t0\t-
            chr1\t1004\t1006\tB\t0\t-
            chr1\t1006\t1009\tC\t0\t-
            """,
            data_size=8,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual([0, 0, 1, 1, 2, 2, 2], as_list(track(Interval("chr1", "+", 1000, 1007, track.refg))))
            self.assertEqual([2, 2, 2, 1, 1, 0, 0], as_list(track(Interval("chr1", "-", 1002, 1009, track.refg))))

        # Valid BED (resolution=2, stranded, scores)
        try_parse(
            2,
            "strand_unaware",
            None,
            """
            chr1\t1000\t1002\tA\t5\t+
            chr1\t1002\t1004\tB\t6\t+
            chr1\t1004\t1007\tC\t7\t+
            chr1\t1002\t1004\tA\t-5\t-
            chr1\t1004\t1006\tB\t-6\t-
            chr1\t1006\t1009\tC\t-7\t-
            """,
            data_size=8,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual([5, 5, 6, 6, 7, 7, 7], as_list(track(Interval("chr1", "+", 1000, 1007, track.refg))))
            self.assertEqual([-7, -7, -7, -6, -6, -5, -5], as_list(track(Interval("chr1", "-", 1002, 1009, track.refg))))

        # Valid BED (resolution=2, mask only -- neither categories nor scores)
        try_parse(
            2,
            "single_stranded",
            None,
            """
            chr1\t1000\t1002
            chr1\t1004\t1007
            """,
            etype='m0',
            data_size=0,
            index_size=2)
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual([1, 1, 0, 0, 1, 1, 1], as_list(track(Interval("chr1", "+", 1000, 1007, track.refg))))

        # m0 mask encoding with name & score columns
        refg = try_parse(1, "single_stranded", None, """chr1\t1000\t1002\tname\t5""", etype='m0')
        self.check_track(Interval('chr1', '+', 999, 1003, refg), [[False], [True], [True], [False]])

        # Unknown chrom
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", ["A", "B", "C"], """chrQ\t1000\t1002\tA\t0\t.""")

        # Unstranded, but '+' or '-' as strand
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", ["A", "B", "C"], """chr1\t1000\t1005\tA\t0\t+""")

        # Stranded, but with '.' as strand
        with self.assertRaises(ValueError):
            try_parse(1, "strand_unaware", ["A", "B", "C"], """chr1\t1000\t1005\tA\t0\t.""")

        # Stranded, but no strand column
        with self.assertRaises(ValueError):
            try_parse(1, "strand_unaware", ["A", "B", "C"], """chr1\t1000\t1005\tA\t0""")

        # Too few columns
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", None, """chr1\t0""")

        # Non-numeric start/end
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", None, """chr1\t@0\t10""")

        # Negative start
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", None, """chr1\t-50\t10""")

        # End <= start
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", None, """chr1\t10\t9""")
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", None, """chr1\t10\t10""")

        # Misaligned start
        with self.assertRaises(ValueError):
            try_parse(10, "single_stranded", None, """chr1\t1\t10""")

        # Misaligned end when not last interval in chromosome
        with self.assertRaises(ValueError):
            try_parse(10, "single_stranded", None, """
                chr1\t0\t9
                chr1\t10\t20
                """)

        # value not encodable as f16
        with self.assertRaises(ValueError):
            try_parse(1, "single_stranded", None, """chr1\t1000\t1002\tname\t100000""", etype='f16')

        # f32 track
        refg = try_parse(1, "single_stranded", None, """chr1\t1000\t1002\tname\t100000""", etype='f32')
        self.check_track(Interval('chr1', '+', 999, 1003, refg), [[0.0], [100000.], [100000.], [0.0]])

        # f8 track with clamping
        maxval = 40000
        d = np.linspace(1, maxval, 256).astype(np.float16)
        refg = try_parse(
            1,
            "strand_unaware", None,
            """
            chr1\t565661\t565662\tchr1:565661..565662,+\t1000\t+
            chr1\t565662\t565663\tchr1:565662..565663,+\t2000\t+
            chr1\t565663\t565664\tchr1:565663..565664,+\t5000\t+
            chr1\t565664\t565665\tchr1:565664..565665,+\t10000\t+
            chr1\t565665\t565666\tchr1:565665..565666,+\t20000\t+
            chr1\t565666\t565667\tchr1:565666..565667,+\t40000\t+
            chr1\t565667\t565668\tchr1:565667..565668,+\t50000\t+
            chr1\t565668\t565669\tchr1:565668..565669,+\t70000\t+
            """,
            etype="f8",
            clamping=True,
            dictionary=d)
        fmaxval = float(maxval)
        self.check_track(
            Interval('chr1', '+', 565661, 565669, refg),
            [[d[6]], [d[13]], [d[32]], [d[64]], [d[128]], [fmaxval], [fmaxval], [fmaxval]],
            # [[  942.], [ 2040.], [ 5020.], [10040.], [20080.], [40000.], [40000.], [40000.]]
            dtype=np.float16,
        )

    def test_restriction(self):
        def try_parse(*args, **kwargs):
            return self._try_parse_wig(*args, **kwargs)

        genome = MiniGenome("hg19")
        # Restriction with dim=2 and resolution=2
        try_parse(
            2,
            2,
            """
            fixedStep chrom=chr1 start=1001 step=2 span=2
            1\t2
            3\t4
            5\t6
            7\t8
            9\t10
            fixedStep chrom=chr2 start=1001 step=2 span=2
            -1\t-2
            -3\t-4
            -5\t-6
            -7\t-8
            -9\t-10
            """,
            data_size=3,
            index_size=1,
            default_value=0,
            restriction=Interval("chr2", "-", 1003, 1007, genome))
        with GenomeTrack(self.tmpfile) as track:
            # Check empty chromosome 1
            self.assertTrue(np.array_equal(np.zeros((10, 2)), track(Interval("chr1", "+", 1000, 1010, track.refg))))

            # Check chromosome 2 and make sure that the data interval
            # got 2 bp cropped off of each end by the restriction.
            a = -np.repeat(np.arange(1, 11).reshape((5, 2)), 2, axis=0)
            ref = np.zeros((10, 2))
            ref[2:8] = a[2:8]
            self.assertTrue(np.array_equal(ref, track(Interval("chr2", "+", 1000, 1010, genome))))

        try_parse(
            2,
            2,
            """
            fixedStep chrom=chr1 start=1 step=2 span=2
            1\t2
            3\t4
            5\t6
            7\t8
            9\t10
            fixedStep chrom=chr1 start=1001 step=2 span=2
            1\t2
            3\t4
            5\t6
            7\t8
            9\t10
            fixedStep chrom=chr1 start=2001 step=2 span=2
            1\t2
            3\t4
            5\t6
            7\t8
            9\t10
            fixedStep chrom=chr1 start=3001 step=2 span=2
            1\t2
            3\t4
            5\t6
            7\t8
            9\t10
            fixedStep chrom=chr1 start=4001 step=2 span=2
            1\t2
            3\t4
            5\t6
            7\t8
            9\t10
            """,
            data_size=7,
            index_size=3,
            default_value=0,
            restriction=Interval("chr1", "-", 1009, 3001, genome))
        with GenomeTrack(self.tmpfile) as track:
            # Check that the following occurred:
            #  1. The 1st interval was excluded entirely.
            #  2. The 2nd interval was cropped away except for the last datum.
            #  3. The 3rd interval was retained entirely
            #  4. The 4th interval was cropped away except for the first datum.
            #  5. The 5th interval was excluded entirely.
            a = np.repeat(np.arange(1, 11).reshape((5, 2)), 2, axis=0)
            ref = np.zeros((4010, 2))
            ref[1008:1010] = a[8:10]
            ref[2000:2010] = a
            ref[3000:3002] = a[0:2]
            self.assertTrue(np.array_equal(ref, track(Interval("chr1", "+", 0, 4010, genome))))

    def test_overlaps(self):
        builder = self.make_builder()

        overlapped = Interval('chr1', '+', 30, 40, builder.refg)
        far_after = Interval('chr1', '+', 60, 70, builder.refg)
        just_before = Interval('chr1', '+', 20, 30, builder.refg)
        just_after = Interval('chr1', '+', 40, 50, builder.refg)
        inside_before = Interval('chr1', '+', 20, 31, builder.refg)
        inside_after = Interval('chr1', '+', 39, 50, builder.refg)

        builder.set_data(inside_before, np.ones(len(inside_before), np.float16))
        with self.assertRaises(ValueError):
            builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(inside_after, np.ones(len(inside_after), np.float16))
        with self.assertRaises(ValueError):
            builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(just_after, np.ones(len(just_after), np.float16))
        builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(just_before, np.ones(len(just_before), np.float16))
        builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(inside_before, np.ones(len(inside_before), np.float16))
        builder.set_data(inside_after, np.ones(len(inside_after), np.float16))
        with self.assertRaises(ValueError):
            builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(just_after, np.ones(len(just_after), np.float16))
        builder.set_data(inside_before, np.ones(len(inside_before), np.float16))
        with self.assertRaises(ValueError):
            builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(far_after, np.ones(len(far_after), np.float16))
        builder.set_data(inside_after, np.ones(len(inside_after), np.float16))
        with self.assertRaises(ValueError):
            builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

        builder = self.make_builder()
        builder.set_data(just_after, np.ones(len(just_after), np.float16))
        builder.set_data(just_before, np.ones(len(just_before), np.float16))
        builder.set_data(overlapped, np.ones(len(overlapped), np.float16))
        del builder

    def test_sense_strand(self):
        builder = self.make_builder(strandedness="strand_unaware")

        interval = Interval('chr1', '-', 10, 15, builder.refg)
        d = np.arange(1, 6, dtype=np.float16)
        self.assertEqual(as_list(d), [1., 2., 3., 4., 5.])

        builder.set_data(interval, d)
        builder.finalize()
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual(as_list(track(interval)), [5., 4., 3., 2., 1.])

        builder = self.make_builder(strandedness="strand_aware")
        builder.set_data(interval, d)
        builder.finalize()
        with GenomeTrack(self.tmpfile) as track:
            self.assertEqual(as_list(track(interval)), [1., 2., 3., 4., 5.])

class TestHistogram(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # intervals to be spanned by data
        # mini genomes chromosome sizes are all 0 except chr2:74682101-74692600
        genome = MiniGenome('hg19')
        cls.intervals_unstranded = [
            Interval('chr2', '+', 0, 5, genome),  # Start of MiniGenome chr2
            Interval('chr2', '+', 5000, 5003, genome),  # Middle of MiniGenome chr2
            Interval('chr2', '+', 10495, 10500, genome),  # End of MiniGenome chr2
        ]
        cls.intervals_stranded = cls.intervals_unstranded + [
            Interval.as_negative_strand(x) for x in cls.intervals_unstranded
        ]

        # 2-dimensional track data for each interval
        cls.data_unstranded = [
            np.array(
                [
                    [0, 0],  # For intervals[0]; -0.0 should get treated as +0.0 in histogram
                    [0, 1],
                    [1, 1],
                    [1, np.inf],
                    [1, np.nan]
                ],
                np.float16),  # this nan should contribute to the nan count
            np.array(
                [
                    [3, 2],  # For intervals[1]
                    [3, 4],
                    [4, 4]
                ],
                np.float16),
            np.array(
                [
                    [4, 4],  # For intervals[2]
                    [4, 5],
                    [5, 6],
                    [5, 6],
                    [5, 6]
                ],
                np.float16),
        ]
        cls.data_stranded = cls.data_unstranded + [np.negative(x) for x in cls.data_unstranded]

        # TODO: check transpose in GenomeTrackBuilder, earlier in file

        def build(is_single_stranded, intervals, data):
            # create tempfile
            fd, tmpfile = mkstemp()
            os.close(fd)

            # build unstranded track into tempfile
            builder = make_builder(tmpfile, strandedness="single_stranded" if is_single_stranded else "strand_unaware", reference_genome=genome, dim=2)
            for interval, datum in zip(intervals, data):
                builder.set_data(interval, datum)
            builder.finalize()

            return tmpfile

        # create tempfiles for the stranded and unstranded tracks
        cls.trackfile_unstranded = build(True, cls.intervals_unstranded, cls.data_unstranded)
        cls.trackfile_stranded = build(False, cls.intervals_stranded, cls.data_stranded)

    @classmethod
    def tearDownClass(cls):
        os.unlink(cls.trackfile_unstranded)
        os.unlink(cls.trackfile_stranded)

    def _histogram(self, array, separate_dims=False):
        """Returns a histogram or list of histograms depending on separate_dims.
        Acts as a reference result for the `GenomeTrack.histogram` method,
        where the data is provided as an array (or list of arrays) directly
        instead of as an interval (or list of intervals).
        """
        if isinstance(array, list):
            array = np.vstack(array)

        # Convert numpy array into a dict of { value : count }
        as_dict = lambda data: dict(Counter(data[~np.isnan(data)]))

        if separate_dims:
            return [as_dict(col) for col in array.T]
        else:
            return as_dict(array.ravel())

    def test_histogram(self):
        def check(trackfile, intervals, data):
            track = GenomeTrack(trackfile)

            # single interval
            self.assertEqual(track.histogram(intervals[0], True), self._histogram(data[0], True))
            self.assertEqual(track.histogram(intervals[0], False), self._histogram(data[0], False))

            # list of intervals
            self.assertEqual(track.histogram(intervals, True), self._histogram(data, True))
            self.assertEqual(track.histogram(intervals, False), self._histogram(data, False))

            # genome-wide; same as using all data intervals
            region = MiniGenome(track.refg)
            self.assertEqual(track.histogram(region, True), self._histogram(data, True))
            self.assertEqual(track.histogram(region, False), self._histogram(data, False))

        # test unstranded and stranded tracks modes
        check(self.trackfile_unstranded, self.intervals_unstranded, self.data_unstranded)
        check(self.trackfile_stranded, self.intervals_stranded, self.data_stranded)

        # test invalid `region` argument types
        with GenomeTrack(self.trackfile_unstranded) as track:
            with self.assertRaises(TypeError):
                track.histogram("foo")
            with self.assertRaises(TypeError):
                track.histogram(0)
            with self.assertRaises(TypeError):
                track.histogram(track)

    def test__intervals_for_genome(self):

        # unstranded track: check it returns exactly one interval per chromosome
        track = GenomeTrack(self.trackfile_unstranded)
        genome = MiniGenome(track.refg)
        intervals = track._intervals_for_genome(genome)
        self.assertEqual(set(map(lambda interval: interval.chromosome, intervals)), set(genome.chromosomes))

        # stranded track: check it returns one interval per strand per chromosome
        track = GenomeTrack(self.trackfile_stranded)
        intervals = track._intervals_for_genome(genome)
        self.assertEqual(
            set(map(lambda interval: (interval.chromosome, interval.strand), intervals)),
            set(product(genome.chromosomes, ('+', '-'))))

        # check error raised if mismatched ref genome
        track = GenomeTrack(self.trackfile_unstranded)
        genome = MiniGenome("hg38.p12")
        with self.assertRaises(ValueError):
            track._intervals_for_genome(genome)

    def test__histogram_counters_for_interval(self):
        track = GenomeTrack(self.trackfile_unstranded)
        counters = track._histogram_counters_for_interval(self.intervals_unstranded[0], 2)

        # it returns a n-dimensional array, each with 2**16 entries (one counter for each possible np.float16 value)
        self.assertIsInstance(counters, np.ndarray)
        self.assertEqual(counters.shape, (2, 2**16))

        self.assertEqual(counters[0].sum(), self.data_unstranded[0].shape[0])
        self.assertEqual(counters[1].sum(), self.data_unstranded[0].shape[0])
