# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
import unittest
import gc
from genome_kit import Interval
from genome_kit import Intron

###############################################################


class NewTest(unittest.TestCase):
    def test(self):
        # Create single interval one at a time.
        # Check initialization, getattr, repr, str, and del
        # Repeat several times to make sure allocation/refcounting works.
        for i in range(10):
            c = Interval("chr5", "+", 2000, 3000, "hg19")
            self.assertEqual(c.chrom, "chr5")
            self.assertEqual(c.chromosome, "chr5")
            self.assertEqual(c.strand, "+")
            self.assertTrue(c.is_positive_strand())
            self.assertEqual(c.start, 2000)
            self.assertEqual(c.end, 3000)
            self.assertEqual(c.end5, Interval("chr5", "+", 2000, 2000, "hg19"))
            self.assertEqual(c.end3, Interval("chr5", "+", 3000, 3000, "hg19"))
            self.assertEqual(c.refg, "hg19")
            self.assertEqual(c.sys, "hg19:chr5")
            self.assertIsInstance(c.__repr__(), str)
            self.assertIsInstance(c.__str__(), str)
            del c  # Make sure single dealloc works
            gc.collect()

        # Create many intervals and check initialization of each after the fact
        tuples = []
        intervals = []
        for i in range(500):
            chrom = "chr{}".format(((i // 7) % 22) + 1)
            strand = ["+", "-"][(i // 5) % 2]
            start = 20000 + i
            end = 20000 + 5 * i
            refg = ["hg19", "hg38.p12"][(i // 3) % 2]
            tuples.append((chrom, strand, start, end, refg))
            intervals.append(Interval(chrom, strand, start, end, refg))

        for (chrom, strand, start, end, refg), c in zip(tuples, intervals):
            self.assertEqual(chrom, c.chrom)
            self.assertEqual(chrom, c.chromosome)
            self.assertEqual(strand, c.strand)
            self.assertEqual(start, c.start)
            self.assertEqual(end, c.end)
            end5 = c.start if strand == "+" else c.end
            end3 = c.end if strand == "+" else c.start
            self.assertEqual(c.end5, Interval(chrom, strand, end5, end5, refg))
            self.assertEqual(c.end3, Interval(chrom, strand, end3, end3, refg))
            self.assertEqual(c.is_positive_strand(), strand == "+")
            if strand == '+':
                self.assertEqual(start, c.end5.start)
                self.assertEqual(end, c.end3.start)
            else:
                self.assertEqual(start, c.end3.start)
                self.assertEqual(end, c.end5.start)
            self.assertEqual(refg, c.refg)
            self.assertEqual(refg, c.reference_genome)
            self.assertEqual(c, Interval(chrom, strand, start, end, refg))
        del tuples
        del intervals
        gc.collect()
        self.assertEqual(gc.garbage, [])


###############################################################


class IterTest(unittest.TestCase):
    def test(self):
        # Positive strand
        c = Interval("chr5", "+", 2000, 2011, "hg19")
        self.assertEqual(len(c), 11)
        self.assertIn(Interval("chr5", "+", 2000, 2001, "hg19"), c)
        self.assertIn(Interval("chr5", "+", 2005, 2006, "hg19"), c)
        self.assertIn(Interval("chr5", "+", 2010, 2011, "hg19"), c)
        self.assertNotIn(Interval("chr5", "+", 2011, 2012, "hg19"), c)
        self.assertNotIn(Interval("chr5", "+", 1999, 2000, "hg19"), c)
        with self.assertRaises(TypeError):
            c[2]  # Does not support indexing
        del c
        gc.collect()
        self.assertEqual(gc.garbage, [])

        # Negative strand
        c = Interval("chr5", "-", 2000, 2011, "hg19")
        self.assertEqual(len(c), 11)
        self.assertIn(Interval("chr5", "-", 2000, 2001, "hg19"), c)
        self.assertIn(Interval("chr5", "-", 2005, 2006, "hg19"), c)
        self.assertIn(Interval("chr5", "-", 2010, 2011, "hg19"), c)
        self.assertNotIn(Interval("chr5", "-", 2011, 2012, "hg19"), c)
        self.assertNotIn(Interval("chr5", "-", 1999, 2000, "hg19"), c)
        with self.assertRaises(TypeError):
            c[2]  # Does not support indexing
        del c
        gc.collect()
        self.assertEqual(gc.garbage, [])


###############################################################


class UniqueTest(unittest.TestCase):
    def test(self):
        intervals = [Interval("chr3", "+", 5000 + i // 2, 8000 + i // 2, "hg19") for i in range(100)]
        unique_intervals = set(intervals)
        self.assertEqual(len(unique_intervals), 50)
        for i in intervals:
            self.assertIn(i, unique_intervals)
        del intervals
        gc.collect()
        self.assertEqual(gc.garbage, [])


###############################################################


class ContainsWithinTest(unittest.TestCase):
    def test(self):
        # yapf: disable
        # The test intervals are in the pattern below
        #   0123456789..
        #       AAAA
        #   BBBB    CCCC
        #    DDDD  EEEE
        #        FF
        #      XY  ZW
        a = Interval("chr1", "+",  4,  8, "hg19")
        b = Interval("chr1", "+",  0,  4, "hg19")
        c = Interval("chr1", "+",  8, 12, "hg19")
        d = Interval("chr1", "+",  1,  5, "hg19")
        e = Interval("chr1", "+",  7, 11, "hg19")
        f = Interval("chr1", "+",  5,  7, "hg19")
        x = Interval("chr1", "+",  3,  4, "hg19")
        y = Interval("chr1", "+",  4,  5, "hg19")
        z = Interval("chr1", "+",  7,  8, "hg19")
        w = Interval("chr1", "+",  8,  9, "hg19")
        items = [a, b, c, d, e, f, x, y, z, w]

        # Filters and version of the filters that operate on opposite strand
        def upstream_of(i):   return lambda j: j.upstream_of(i)
        def upstream_of_r(i): return lambda j: j.as_opposite_strand().upstream_of(i.as_opposite_strand())
        def dnstream_of(i):   return lambda j: j.dnstream_of(i)
        def dnstream_of_r(i): return lambda j: j.as_opposite_strand().dnstream_of(i.as_opposite_strand())
        def contains(i):      return lambda j: j.contains(i)
        def contains_r(i):    return lambda j: j.as_opposite_strand().contains(i.as_opposite_strand())
        def _contains(i):     return lambda j: i in j
        def _contains_r(i):   return lambda j: i.as_opposite_strand() in j.as_opposite_strand()
        def within(i):        return lambda j: j.within(i)
        def within_r(i):      return lambda j: j.as_opposite_strand().within(i.as_opposite_strand())
        def overlaps(i):      return lambda j: j.overlaps(i)
        def overlaps_r(i):    return lambda j: j.as_opposite_strand().overlaps(i.as_opposite_strand())

        self.assertEqual(list(filter(upstream_of(a),   items)), [b, x])
        self.assertEqual(list(filter(dnstream_of_r(a), items)), [b, x])
        self.assertEqual(list(filter(dnstream_of(a),   items)), [c, w])
        self.assertEqual(list(filter(upstream_of_r(a), items)), [c, w])
        self.assertEqual(list(filter(contains(a),    items)), [a])
        self.assertEqual(list(filter(contains_r(a),  items)), [a])
        self.assertEqual(list(filter(_contains(a),   items)), [a])
        self.assertEqual(list(filter(_contains_r(a), items)), [a])
        self.assertEqual(list(filter(within(a),      items)), [a, f, y, z])
        self.assertEqual(list(filter(within_r(a),    items)), [a, f, y, z])
        self.assertEqual(list(filter(overlaps(a),    items)), [a, d, e, f, y, z])
        self.assertEqual(list(filter(overlaps_r(a),  items)), [a, d, e, f, y, z])

        self.assertEqual(list(filter(upstream_of(b),   items)), [])
        self.assertEqual(list(filter(dnstream_of_r(b), items)), [])
        self.assertEqual(list(filter(dnstream_of(b),   items)), [a, c, e, f, y, z, w])
        self.assertEqual(list(filter(upstream_of_r(b), items)), [a, c, e, f, y, z, w])
        self.assertEqual(list(filter(contains(b),    items)), [b])
        self.assertEqual(list(filter(contains_r(b),  items)), [b])
        self.assertEqual(list(filter(_contains(b),   items)), [b])
        self.assertEqual(list(filter(_contains_r(b), items)), [b])
        self.assertEqual(list(filter(within(b),      items)), [b, x])
        self.assertEqual(list(filter(within_r(b),    items)), [b, x])
        self.assertEqual(list(filter(overlaps(b),    items)), [b, d, x])
        self.assertEqual(list(filter(overlaps_r(b),  items)), [b, d, x])

        self.assertEqual(list(filter(upstream_of(c),   items)), [a, b, d, f, x, y, z])
        self.assertEqual(list(filter(dnstream_of_r(c), items)), [a, b, d, f, x, y, z])
        self.assertEqual(list(filter(dnstream_of(c),   items)), [])
        self.assertEqual(list(filter(upstream_of_r(c), items)), [])
        self.assertEqual(list(filter(contains(c),    items)), [c])
        self.assertEqual(list(filter(contains_r(c),  items)), [c])
        self.assertEqual(list(filter(_contains(c),   items)), [c])
        self.assertEqual(list(filter(_contains_r(c), items)), [c])
        self.assertEqual(list(filter(within(c),      items)), [c, w])
        self.assertEqual(list(filter(within_r(c),    items)), [c, w])
        self.assertEqual(list(filter(overlaps(c),    items)), [c, e, w])
        self.assertEqual(list(filter(overlaps_r(c),  items)), [c, e, w])

        self.assertEqual(list(filter(upstream_of(d),   items)), [])
        self.assertEqual(list(filter(dnstream_of_r(d), items)), [])
        self.assertEqual(list(filter(dnstream_of(d),   items)), [c, e, f, z, w])
        self.assertEqual(list(filter(upstream_of_r(d), items)), [c, e, f, z, w])
        self.assertEqual(list(filter(contains(d),    items)), [d])
        self.assertEqual(list(filter(contains_r(d),  items)), [d])
        self.assertEqual(list(filter(_contains(d),   items)), [d])
        self.assertEqual(list(filter(_contains_r(d), items)), [d])
        self.assertEqual(list(filter(within(d),      items)), [d, x, y])
        self.assertEqual(list(filter(within_r(d),    items)), [d, x, y])
        self.assertEqual(list(filter(overlaps(d),    items)), [a, b, d, x, y])
        self.assertEqual(list(filter(overlaps_r(d),  items)), [a, b, d, x, y])

        self.assertEqual(list(filter(upstream_of(e),   items)), [b, d, f, x, y])
        self.assertEqual(list(filter(dnstream_of_r(e), items)), [b, d, f, x, y])
        self.assertEqual(list(filter(dnstream_of(e),   items)), [])
        self.assertEqual(list(filter(upstream_of_r(e), items)), [])
        self.assertEqual(list(filter(contains(e),    items)), [e])
        self.assertEqual(list(filter(contains_r(e),  items)), [e])
        self.assertEqual(list(filter(_contains(e),   items)), [e])
        self.assertEqual(list(filter(_contains_r(e), items)), [e])
        self.assertEqual(list(filter(within(e),      items)), [e, z, w])
        self.assertEqual(list(filter(within_r(e),    items)), [e, z, w])
        self.assertEqual(list(filter(overlaps(e),    items)), [a, c, e, z, w])
        self.assertEqual(list(filter(overlaps_r(e),  items)), [a, c, e, z, w])

        self.assertEqual(list(filter(upstream_of(f),   items)), [b, d, x, y])
        self.assertEqual(list(filter(dnstream_of_r(f), items)), [b, d, x, y])
        self.assertEqual(list(filter(dnstream_of(f),   items)), [c, e, z, w])
        self.assertEqual(list(filter(upstream_of_r(f), items)), [c, e, z, w])
        self.assertEqual(list(filter(contains(f),    items)), [a, f])
        self.assertEqual(list(filter(contains_r(f),  items)), [a, f])
        self.assertEqual(list(filter(_contains(f),   items)), [a, f])
        self.assertEqual(list(filter(_contains_r(f), items)), [a, f])
        self.assertEqual(list(filter(within(f),      items)), [f])
        self.assertEqual(list(filter(within_r(f),    items)), [f])
        self.assertEqual(list(filter(overlaps(f),    items)), [a, f])
        self.assertEqual(list(filter(overlaps_r(f),  items)), [a, f])

        self.assertEqual(list(filter(upstream_of(x),   items)), [])
        self.assertEqual(list(filter(dnstream_of_r(x), items)), [])
        self.assertEqual(list(filter(dnstream_of(x),   items)), [a, c, e, f, y, z, w])
        self.assertEqual(list(filter(upstream_of_r(x), items)), [a, c, e, f, y, z, w])
        self.assertEqual(list(filter(contains(x),    items)), [b, d, x])
        self.assertEqual(list(filter(contains_r(x),  items)), [b, d, x])
        self.assertEqual(list(filter(_contains(x),   items)), [b, d, x])
        self.assertEqual(list(filter(_contains_r(x), items)), [b, d, x])
        self.assertEqual(list(filter(within(x),      items)), [x])
        self.assertEqual(list(filter(within_r(x),    items)), [x])
        self.assertEqual(list(filter(overlaps(x),    items)), [b, d, x])
        self.assertEqual(list(filter(overlaps_r(x),  items)), [b, d, x])

        self.assertEqual(list(filter(upstream_of(y),   items)), [b, x])
        self.assertEqual(list(filter(dnstream_of_r(y), items)), [b, x])
        self.assertEqual(list(filter(dnstream_of(y),   items)), [c, e, f, z, w])
        self.assertEqual(list(filter(upstream_of_r(y), items)), [c, e, f, z, w])
        self.assertEqual(list(filter(contains(y),    items)), [a, d, y])
        self.assertEqual(list(filter(contains_r(y),  items)), [a, d, y])
        self.assertEqual(list(filter(_contains(y),   items)), [a, d, y])
        self.assertEqual(list(filter(_contains_r(y), items)), [a, d, y])
        self.assertEqual(list(filter(within(y),      items)), [y])
        self.assertEqual(list(filter(within_r(y),    items)), [y])
        self.assertEqual(list(filter(overlaps(y),    items)), [a, d, y])
        self.assertEqual(list(filter(overlaps_r(y),  items)), [a, d, y])

        self.assertEqual(list(filter(upstream_of(z),   items)), [b, d, f, x, y])
        self.assertEqual(list(filter(dnstream_of_r(z), items)), [b, d, f, x, y])
        self.assertEqual(list(filter(dnstream_of(z),   items)), [c, w])
        self.assertEqual(list(filter(upstream_of_r(z), items)), [c, w])
        self.assertEqual(list(filter(contains(z),    items)), [a, e, z])
        self.assertEqual(list(filter(contains_r(z),  items)), [a, e, z])
        self.assertEqual(list(filter(_contains(z),   items)), [a, e, z])
        self.assertEqual(list(filter(_contains_r(z), items)), [a, e, z])
        self.assertEqual(list(filter(within(z),      items)), [z])
        self.assertEqual(list(filter(within_r(z),    items)), [z])
        self.assertEqual(list(filter(overlaps(z),    items)), [a, e, z])
        self.assertEqual(list(filter(overlaps_r(z),  items)), [a, e, z])

        self.assertEqual(list(filter(upstream_of(w),   items)), [a, b, d, f, x, y, z])
        self.assertEqual(list(filter(dnstream_of_r(w), items)), [a, b, d, f, x, y, z])
        self.assertEqual(list(filter(dnstream_of(w),   items)), [])
        self.assertEqual(list(filter(upstream_of_r(w), items)), [])
        self.assertEqual(list(filter(contains(w),    items)), [c, e, w])
        self.assertEqual(list(filter(contains_r(w),  items)), [c, e, w])
        self.assertEqual(list(filter(_contains(w),   items)), [c, e, w])
        self.assertEqual(list(filter(_contains_r(w), items)), [c, e, w])
        self.assertEqual(list(filter(within(w),      items)), [w])
        self.assertEqual(list(filter(within_r(w),    items)), [w])
        self.assertEqual(list(filter(overlaps(w),    items)), [c, e, w])
        self.assertEqual(list(filter(overlaps_r(w),  items)), [c, e, w])
        # yapf: enable


########################################################


class TestInterval(unittest.TestCase):
    def test_initializers(self):
        # Initialize with string identifier for strand
        self.assertIsInstance(Interval('chr1', '+', 9, 9, 'hg19'), Interval)
        self.assertIsInstance(Interval('chr1', '-', 9, 9, 'hg19'), Interval)
        self.assertIsInstance(Interval('chr1', '+', 9, 20, 'hg19'), Interval)
        self.assertIsInstance(Interval('chr1', '-', 9, 20, 'hg19'), Interval)
        self.assertIsInstance(Interval.from_dna0('chr1', '-', 9, 20, 'hg19'), Interval)

        with self.assertRaises(TypeError):
            Interval('chr1', 3, 9, 20, 'hg19')

        # Test invalid start and end values
        with self.assertRaises(ValueError):
            Interval('chr1', '+', 10, 5, 'hg19')

        with self.assertRaises(ValueError):
            Interval('chr99', '+', 5, 10, 'hg19')

        with self.assertRaises(ValueError):
            Interval.from_rna1('chr1', -100, -200, 'hg19')

        with self.assertRaises(ValueError):
            Interval.from_rna1('chr1', -100, 200, 'hg19')

        self.assertIsInstance(str(Interval('chr1', '+', 100, 200, 'hg19')), str)

    def test_forward_strand(self):
        from_dna0 = Interval.from_dna0('chr1', '+', 9, 20, 'hg19')
        from_dna1 = Interval.from_dna1('chr1', '+', 10, 20, 'hg19')
        from_rna1 = Interval.from_rna1('chr1', 10, 20, 'hg19')
        derived_from_dna0 = Intron.from_dna0('chr1', '+', 9, 20, 'hg19')
        derived_from_dna1 = Intron.from_dna1('chr1', '+', 10, 20, 'hg19')
        derived_from_rna1 = Intron.from_rna1('chr1', 10, 20, 'hg19')

        self.assertEqual(len(from_dna0), 11)

        self.assertEqual(from_dna0.as_dna0(), (9, 20))
        self.assertEqual(from_dna0.as_dna1(), (10, 20))
        self.assertEqual(from_dna0.as_rna1(), (10, 20))
        self.assertEqual(from_dna0.as_ucsc(), 'chr1:10-20')

        self.assertEqual(from_dna0, from_dna1)
        self.assertEqual(from_dna0, from_rna1)

        self.assertEqual(from_dna0, derived_from_dna0)
        self.assertEqual(from_dna1, derived_from_dna1)
        self.assertEqual(from_rna1, derived_from_rna1)

    def test_reverse_strand(self):
        from_dna0 = Interval('chr1', '-', 9, 20, 'hg19')
        from_dna1 = Interval.from_dna1('chr1', '-', 10, 20, 'hg19')
        from_rna1 = Interval.from_rna1('chr1', -20, -10, 'hg19')

        self.assertEqual(len(from_dna0), 11)

        self.assertEqual(from_dna0.as_dna0(), (9, 20))
        self.assertEqual(from_dna0.as_dna1(), (10, 20))
        self.assertEqual(from_dna0.as_rna1(), (-20, -10))
        self.assertEqual(from_dna0.as_ucsc(), 'chr1:10-20')

        self.assertEqual(from_dna0, from_dna1)
        self.assertEqual(from_dna0, from_rna1)

    def test_from_coord(self):
        # Positive strand
        coord = Interval.from_dna0_coord('chr1', '+', 15, 'hg19')
        derived_coord = Intron.from_dna0_coord('chr1', '+', 15, 'hg19')
        self.assertEqual(coord.start, 15)
        self.assertEqual(coord.end, 16)
        self.assertEqual(coord, derived_coord)

        # Negative strand
        coord = Interval.from_dna0_coord('chr1', '-', 15, 'hg19')
        self.assertEqual(coord.start, 15)
        self.assertEqual(coord.end, 16)

    def test_interval_attr(self):
        interval = Interval('chr1', '+', 9, 9, 'hg19')
        self.assertIs(interval.interval, interval)

    def test_shifted_by(self):
        gap = Interval('chr1', '+', 15, 15, 'hg19', 30, 31)

        interval = gap.shift(5)
        self.assertEqual(interval.chrom, gap.chrom)
        self.assertEqual(interval.strand, gap.strand)
        self.assertEqual(interval.refg, gap.refg)
        self.assertEqual(interval.start, 20)
        self.assertEqual(interval.end, 20)
        self.assertEqual(interval.anchor, 30)
        self.assertEqual(interval.anchor_offset, 31)

        interval = gap.shift(-5)
        self.assertEqual(interval.start, 10)
        self.assertEqual(interval.end, 10)
        self.assertEqual(interval.anchor, 30)
        self.assertEqual(interval.anchor_offset, 31)

        interval = gap.as_negative_strand().shift(5)
        self.assertEqual(interval.start, 10)
        self.assertEqual(interval.end, 10)
        self.assertEqual(interval.anchor, 30)
        self.assertEqual(interval.anchor_offset, 31)

    def test_expanded_by(self):
        gap = Interval('chr1', '+', 15, 15, 'hg19', 30, 31)

        interval = gap.expand(5)
        self.assertEqual(interval.chrom, gap.chrom)
        self.assertEqual(interval.strand, gap.strand)
        self.assertEqual(interval.refg, gap.refg)
        self.assertEqual(interval.start, 10)
        self.assertEqual(interval.end, 20)
        self.assertEqual(interval.anchor, 30)
        self.assertEqual(interval.anchor_offset, 31)

        interval = gap.expand(10, 5)
        self.assertEqual(interval.start, 5)
        self.assertEqual(interval.end, 20)
        self.assertEqual(interval.anchor, 30)
        self.assertEqual(interval.anchor_offset, 31)

        interval = gap.as_negative_strand().expand(10, 5)
        self.assertEqual(interval.start, 10)
        self.assertEqual(interval.end, 25)
        self.assertEqual(interval.anchor, 30)
        self.assertEqual(interval.anchor_offset, 31)

        with self.assertRaises(ValueError):
            gap.expand(-1, 0)

    def test_end5_end3(self):
        interval = Interval('chr1', '+', 10, 15, 'hg19', 30, 31)
        self.assertEqual(interval.end5, Interval('chr1', '+', 10, 10, 'hg19', 30, 31))
        self.assertEqual(interval.end3, Interval('chr1', '+', 15, 15, 'hg19', 30, 31))
        self.assertEqual(interval.end5.anchor, 30)
        self.assertEqual(interval.end3.anchor_offset, 31)

        interval = Interval('chr1', '-', 10, 15, 'hg19', 30, 31)
        self.assertEqual(interval.end5, Interval('chr1', '-', 15, 15, 'hg19', 30, 31))
        self.assertEqual(interval.end3, Interval('chr1', '-', 10, 10, 'hg19', 30, 31))
        self.assertEqual(interval.end5.anchor, 30)
        self.assertEqual(interval.end3.anchor_offset, 31)

    def test_comparisons_forward(self):
        a = Interval.from_rna1('chr1', 10, 20, 'hg19')
        b = Interval.from_rna1('chr1', 5, 15, 'hg19')
        c = Interval.from_rna1('chr1', 5, 8, 'hg19')
        d = Interval.from_rna1('chr1', 12, 17, 'hg19')

        self.assertTrue(c.upstream_of(a))
        self.assertFalse(b.upstream_of(a))

        self.assertTrue(a.dnstream_of(c))
        self.assertFalse(a.dnstream_of(b))

        self.assertTrue(a.contains(d))
        self.assertFalse(a.contains(b))

        self.assertTrue(a.upstream_of(Interval.from_rna1('chr1', 21, 21, 'hg19')))
        self.assertFalse(a.upstream_of(Interval.from_rna1('chr1', 20, 20, 'hg19')))
        self.assertFalse(a.upstream_of(Interval.from_rna1('chr1', 19, 19, 'hg19')))

        self.assertTrue(a.contains(Interval.from_rna1('chr1', 15, 15, 'hg19')))

        self.assertTrue(a.dnstream_of(Interval.from_rna1('chr1', 9, 9, 'hg19')))
        self.assertFalse(a.dnstream_of(Interval.from_rna1('chr1', 10, 10, 'hg19')))

        # Test incompatible reference genomes
        e = Interval.from_rna1('chr1', 12, 17, 'hg38.p12')
        with self.assertRaises(ValueError):
            d.upstream_of(e)
        with self.assertRaises(ValueError):
            d.dnstream_of(e)
        with self.assertRaises(ValueError):
            d.contains(e)
        with self.assertRaises(ValueError):
            d.within(e)
        with self.assertRaises(ValueError):
            d.overlaps(e)
        with self.assertRaises(ValueError):
            d == e
        with self.assertRaises(ValueError):
            d != e

    def test_comparisons_reverse(self):
        a = Interval.from_rna1('chr1', -20, -10, 'hg19')
        b = Interval.from_rna1('chr1', -15, -5, 'hg19')
        c = Interval.from_rna1('chr1', -8, -5, 'hg19')
        d = Interval.from_rna1('chr1', -17, -12, 'hg19')

        self.assertTrue(a.upstream_of(c))
        self.assertFalse(a.upstream_of(b))
        self.assertTrue(a.overlaps(b))

        self.assertTrue(c.dnstream_of(a))
        self.assertFalse(b.dnstream_of(a))

        self.assertTrue(a.contains(d))
        self.assertFalse(a.contains(b))

        self.assertTrue(d.within(a))
        self.assertFalse(b.within(a))

        self.assertTrue(a.upstream_of(Interval.from_rna1('chr1', -9, -9, 'hg19')))
        self.assertFalse(a.upstream_of(Interval.from_rna1('chr1', -10, -10, 'hg19')))
        self.assertFalse(a.upstream_of(Interval.from_rna1('chr1', -11, -11, 'hg19')))

        self.assertTrue(a.dnstream_of(Interval.from_rna1('chr1', -21, -21, 'hg19')))
        self.assertFalse(a.dnstream_of(Interval.from_rna1('chr1', -20, -20, 'hg19')))

    def test_input_validation(self):

        with self.assertRaises(ValueError):
            Interval('chr1', '+', 1, -5, 'hg19')

        with self.assertRaises(ValueError):
            Interval('chr1', '+', 10, 5, 'hg19')

        with self.assertRaises(TypeError):
            Interval('chr1', 2, 1, 10, 'hg19')

        with self.assertRaises(ValueError):
            Interval('chr99', '+', -1, 1, 'hg19')

        with self.assertRaises(ValueError):
            Interval.from_rna1('chr1', 10, 5, 'hg19')

        with self.assertRaises(ValueError):
            Interval.from_rna1('chr1', -5, 5, 'hg19')

    def test_change_strand(self):
        a = Interval('chr1', '+', 5, 10, 'hg19')
        b = a.as_opposite_strand()
        c = b.as_opposite_strand()

        self.assertEqual(a, c)
        self.assertEqual(a, Interval('chr1', '+', 5, 10, 'hg19'))
        self.assertEqual(b, Interval('chr1', '-', 5, 10, 'hg19'))
        self.assertEqual(a, a.as_positive_strand())
        self.assertNotEqual(a, a.as_negative_strand())
        self.assertNotEqual(b, b.as_positive_strand())
        self.assertEqual(b, b.as_negative_strand())

    def test_anchor(self):
        interval = Interval('chr1', '+', 5, 10, 'hg19', 'start')
        self.assertEqual(interval.anchor, 5)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '+', 5, 10, 'hg19', 'end')
        self.assertEqual(interval.anchor, 10)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '+', 5, 10, 'hg19', '5p')
        self.assertEqual(interval.anchor, 5)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '-', 5, 10, 'hg19', '5p')
        self.assertEqual(interval.anchor, 10)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '+', 5, 10, 'hg19', '3p')
        self.assertEqual(interval.anchor, 10)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '-', 5, 10, 'hg19', '3p')
        self.assertEqual(interval.anchor, 5)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '+', 5, 10, 'hg19', 'center')
        self.assertEqual(interval.anchor, 7)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '-', 5, 10, 'hg19', 'center')
        self.assertEqual(interval.anchor, 8)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '+', 5, 10, 'hg19', 9)
        self.assertEqual(interval.anchor, 9)
        self.assertEqual(interval.anchor_offset, 0)

        interval = Interval('chr1', '+', 5, 10, 'hg19', 9, 3)
        self.assertEqual(interval.anchor, 9)
        self.assertEqual(interval.anchor_offset, 3)

        interval = Interval('chr1', '+', 5, 10, 'hg19', Interval('chr1', '+', 10, 10, 'hg19'), 7)
        self.assertEqual(interval.anchor, 10)
        self.assertEqual(interval.anchor_offset, 7)

        interval = Interval('chr1', '+', 5, 10, 'hg19', 5, 2)
        anchored = interval.with_anchor(20, 22)
        self.assertEqual(anchored.chrom, interval.chrom)
        self.assertEqual(anchored.start, interval.start)
        self.assertEqual(anchored.end, interval.end)
        self.assertEqual(anchored.refg, interval.refg)
        self.assertEqual(anchored.anchor, 20)
        self.assertEqual(anchored.anchor_offset, 22)

        with self.assertRaises(ValueError):  # Nonsense anchor string
            Interval('chr1', '+', 5, 10, 'hg19', 'foo')

        with self.assertRaises(ValueError):  # Different reference genome
            Interval('chr1', '+', 5, 10, 'hg19', Interval('chr1', '+', 10, 10, 'hg38.p12'))

        with self.assertRaises(ValueError):  # Different chromosome
            Interval('chr1', '+', 5, 10, 'hg19', Interval('chr2', '+', 10, 10, 'hg19'))

        with self.assertRaises(ValueError):  # Different strand
            Interval('chr1', '+', 5, 10, 'hg19', Interval('chr1', '-', 10, 10, 'hg19'))

    def test_serialize(self):
        intervals = [
            Interval('chr1', '+', 10, 20, 'hg19'),
            Interval('chr8', '+', 100, 200, 'hg38.p12'),
            Interval('chrX', '-', 10, 200, 'hg19'),
            Interval('chr2', '-', 10, 20, 'hg19', 20),
            Interval('chrY', '-', 10, 20, 'hg19', 20, 25),
        ]

        # Check that pickling / unpickling works
        import pickle
        stored = pickle.dumps(intervals, pickle.HIGHEST_PROTOCOL)
        loaded = pickle.loads(stored)
        self.assertEqual(loaded, intervals)

    def test_spanning(self):

        # success case
        interval1 = Interval('chr1', '+', 10, 20, 'hg19')
        interval2 = Interval('chr1', '+', 30, 40, 'hg19')

        spanning = Interval.spanning(interval1, interval2)
        derived_spanning = Intron.spanning(interval1, interval2)
        self.assertEqual('chr1', spanning.chromosome)
        self.assertEqual('+', spanning.strand)
        self.assertEqual(10, spanning.start)
        self.assertEqual(40, spanning.end)
        self.assertEqual(interval1.reference_genome, spanning.reference_genome)
        self.assertEqual(spanning, derived_spanning)

        # error cases, different chromosomes
        interval3 = Interval('chr2', '+', 30, 40, 'hg19')
        with self.assertRaises(ValueError):
            Interval.spanning(interval1, interval3)

        # different strand
        interval4 = Interval('chr1', '-', 30, 40, 'hg19')
        with self.assertWarns(UserWarning):
            Interval.spanning(interval1, interval4)

        # different reference genome
        interval5 = Interval('chr1', '+', 30, 40, 'hg38.p12')
        with self.assertRaises(ValueError):
            Interval.spanning(interval1, interval5)

        with self.assertRaises(ValueError):
            Interval.spanning(interval1.with_anchor(0), interval2)
        with self.assertRaises(ValueError):
            Interval.spanning(interval1, interval2.with_anchor(0))


    def test_empty_interval(self):
        empty_interval = Interval("chr1", "+", 10, 10, "hg19")

        with self.assertRaises(ValueError):
            empty_interval.as_rna1()

        with self.assertRaises(ValueError):
            empty_interval.as_dna1()

    def test_intersect(self):
        left = Interval('chr1', '+', 10, 20, 'hg19')
        mid = Interval('chr1', '+', 15, 25, 'hg19')
        right = Interval('chr1', '+', 20, 30, 'hg19')
        far_right = Interval('chr1', '+', 50, 55, 'hg19')
        inside = Interval('chr1', '+', 22, 22, 'hg19')

        self.assertEqual(left.intersect(mid), Interval('chr1', '+', 15, 20, 'hg19'))
        self.assertEqual(mid.intersect(left), Interval('chr1', '+', 15, 20, 'hg19'))
        self.assertEqual(right.intersect(mid), Interval('chr1', '+', 20, 25, 'hg19'))
        self.assertEqual(mid.intersect(right), Interval('chr1', '+', 20, 25, 'hg19'))
        self.assertIsNone(left.intersect(far_right))
        self.assertIsNone(left.intersect(right))
        self.assertIsNone(right.intersect(left))
        self.assertEqual(mid.intersect(inside), inside)

        neg_left = left.as_negative_strand()
        neg_mid = mid.as_negative_strand()
        neg_right = right.as_negative_strand()
        neg_far_right = far_right.as_negative_strand()
        neg_inside = inside.as_negative_strand()

        self.assertEqual(neg_left.intersect(neg_mid), Interval('chr1', '-', 15, 20, 'hg19'))
        self.assertEqual(neg_mid.intersect(neg_left), Interval('chr1', '-', 15, 20, 'hg19'))
        self.assertEqual(neg_right.intersect(neg_mid), Interval('chr1', '-', 20, 25, 'hg19'))
        self.assertEqual(neg_mid.intersect(neg_right), Interval('chr1', '-', 20, 25, 'hg19'))
        self.assertIsNone(neg_left.intersect(neg_far_right))
        self.assertIsNone(neg_left.intersect(neg_right))
        self.assertIsNone(neg_right.intersect(neg_left))
        self.assertEqual(neg_mid.intersect(neg_inside), neg_inside)

        self.assertIsNone(left.intersect(Interval(left.chrom, left.strand, left.start, left.end, 'hg38.p12')))
        self.assertIsNone(neg_left.intersect(left))

        with self.assertRaises(ValueError):
            left.intersect(right.with_anchor('5p'))
        with self.assertRaises(ValueError):
            left.with_anchor('5p').intersect(right)

    def test_subtract(self):
        left = Interval('chr1', '+', 10, 20, 'hg19')
        mid = Interval('chr1', '+', 15, 25, 'hg19')
        right = Interval('chr1', '+', 20, 30, 'hg19')
        inside = Interval('chr1', '+', 22, 23, 'hg19')
        covered = Interval('chr1', '+', 10, 30, 'hg19')
        neg_left = left.as_negative_strand()
        diff_refg_left = Interval(left.chrom, left.strand, left.start, left.end, 'hg38.p12')

        self.assertEqual(left.subtract(mid), [Interval('chr1', '+', 10, 15, 'hg19')])
        self.assertEqual(mid.subtract(left), [Interval('chr1', '+', 20, 25, 'hg19')])
        self.assertEqual(right.subtract(mid), [Interval('chr1', '+', 25, 30, 'hg19')])
        self.assertEqual(mid.subtract(right), [Interval('chr1', '+', 15, 20, 'hg19')])
        self.assertEqual(left.subtract(right), [left])
        self.assertEqual(right.subtract(left), [right])
        self.assertEqual(mid.subtract(inside), [Interval('chr1', '+', 15, 22, 'hg19'),
                                                Interval('chr1', '+', 23, 25, 'hg19')])
        self.assertEqual(mid.subtract(covered), [])
        self.assertEqual(left.subtract(neg_left), [left])
        self.assertEqual(left.subtract(diff_refg_left), [left])

        with self.assertRaises(ValueError):
            left.subtract(right.with_anchor('5p'))
        with self.assertRaises(ValueError):
            left.with_anchor('5p').subtract(right)

    def test_midpoint(self):
        point = Interval('chr1', '+', 100, 100, 'hg19')
        self.assertEqual(point.midpoint, point)

        one = point.expand(0, 1)
        self.assertEqual(one.midpoint, point.expand(0, 1))

        two = point.expand(0, 2)
        self.assertEqual(two.midpoint, point.shift(1))

        three = point.expand(0, 3)
        self.assertEqual(three.midpoint, point.shift(1).expand(0, 1))

    def test_distance(self):
        point = Interval('chr1', '+', 100, 100, 'hg19')
        one = point.expand(0, 1)
        two = point.expand(0, 2)
        three = point.expand(0, 3)

        self.assertEqual(point.distance(one), 0.5)
        self.assertEqual(point.distance(two), 1.0)
        self.assertEqual(point.distance(three), 1.5)

        self.assertEqual(point.distance(one, method='end5'), 0.0)
        self.assertEqual(point.distance(two, method='end5'), 0.0)
        self.assertEqual(point.distance(three, method='end5'), 0.0)

        self.assertEqual(point.distance(one, method='end3'), 1.0)
        self.assertEqual(point.distance(two, method='end3'), 2.0)
        self.assertEqual(point.distance(three, method='end3'), 3.0)

        self.assertEqual(point.distance([one, two, three], method='end3'), [1.0, 2.0, 3.0])

        with self.assertRaises(ValueError):
            point.distance(Interval('chr2', point.strand, point.start, point.end, point.refg))
        with self.assertRaises(ValueError):
            point.distance(Interval(point.chrom, point.strand, point.start, point.end, 'hg38.p12'))
        with self.assertRaises(ValueError):
            point.distance(point.with_anchor(0))
        with self.assertWarns(UserWarning):
            point.distance(point.as_opposite_strand())


###############################################################

if __name__ == "__main__":
    unittest.main()
