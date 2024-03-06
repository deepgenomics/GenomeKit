# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
import unittest
from genome_kit import Interval, Variant
from . import MiniGenome


class TestVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        cls.genome = MiniGenome('test_genome')

    @classmethod
    def tearDownClass(cls):
        del cls.genome

    def test_variant_init(self):
        chromosome = 'chr1'
        position = 9
        ref = 'C'
        alt = 'G'

        variant = Variant(chromosome, position, ref, alt, self.genome)

        self.assertEqual(variant.chromosome, chromosome)
        self.assertEqual(variant.position, position)
        self.assertEqual(variant.ref, ref)
        self.assertEqual(variant.alt, alt)
        self.assertEqual(variant.reference_genome, self.genome.reference_genome)

    def test_variant_init_from_string(self):
        variant_str = 'chr1:10:C:G'

        chromosome = 'chr1'
        position = 9
        ref = 'C'
        alt = 'G'

        variant_from_str = Variant.from_string(variant_str, self.genome)
        variant = Variant(chromosome, position, ref, alt, self.genome)

        self.assertEqual(variant_from_str, variant)
        self.assertEqual(variant.as_variant_string(), variant_str)

    def test_ensembl_chromosome_name(self):
        genome37 = MiniGenome('test_genome')
        variant = Variant.from_string("2:1:A:T", genome37)
        self.assertEqual(variant.chromosome, 'chr2')

    def test_preprocess_variants(self):
        variants = u'1:,,,10:c:g'
        variant = Variant.from_string(variants, self.genome)
        self.assertEqual(variant.chromosome, 'chr1')
        self.assertEqual(variant.position, 9)
        self.assertEqual(variant.ref, 'C')
        self.assertEqual(variant.alt, 'G')

    def test_invalid_variant(self):
        with self.assertRaisesRegex(ValueError, 'Invalid variant:'):
            Variant.from_string(4, self.genome)

    def test_too_many_components(self):
        with self.assertRaisesRegex(ValueError, "Variant must have exactly four components:"):
            Variant.from_string('1:2:3:4:5', self.genome)

    def test_invalid_position(self):
        with self.assertRaisesRegex(ValueError, 'Invalid position in variant:'):
            Variant.from_string('chr1:xxx:T:G', self.genome)

    def test_not_match_reference(self):
        with self.assertRaisesRegex(ValueError, 'Variant string was invalid:') as cm:
            Variant.from_string('chr1:2:T:G', self.genome)
        self.assertRegex(str(cm.exception.__cause__), "Variant's ref sequence does not ")

    def test_preprocessing(self):
        variant = Variant.from_string('1:11:-:G', self.genome)
        self.assertEqual(variant.chromosome, 'chr1')
        self.assertEqual(variant.ref, '')

        variant = Variant.from_string('chr1:11:G:.', self.genome)
        self.assertEqual(variant.alt, '')

    def test_normalize_variant(self):
        # Insertion followed by substitution - results in two normalized variants
        variant = Variant('chr1', 16, 'ACGTA', 'TGATACTA', self.genome)
        normalized_variant = variant._normalized_variant
        self.assertEqual(normalized_variant, [Variant('chr1', 16, 'ACGTA', 'TGATA', self.genome), Variant('chr1', 21, '', 'CTA', self.genome)])

        # Insertion followed by substitution - results in two normalized variants
        variant = Variant('chr1', 14, 'GTACGT', 'TTT', self.genome)
        normalized_variant = variant._normalized_variant
        self.assertEqual(normalized_variant, [Variant('chr1', 14, 'GTA', 'TTT', self.genome), Variant('chr1', 17, 'CGT', '', self.genome)])

        # Simple substitution - returns single normalized variant
        variant = Variant('chr1', 9, 'C', 'G', self.genome)
        normalized_variant = variant._normalized_variant
        self.assertEqual(normalized_variant, [Variant('chr1', 9, 'C', 'G', self.genome)])

        # Noop substitution followed by insertion - returns a single normalized insertion variant
        variant = Variant('chr1', 9, 'C', 'CG', self.genome)
        normalized_variant = variant._normalized_variant
        self.assertEqual(normalized_variant, [Variant('chr1', 10, '', 'G', self.genome)])

        # Noop substitution followed by deletion - returns a single normalized deletion variant
        variant = Variant('chr1', 9, 'CG', 'C', self.genome)
        normalized_variant = variant._normalized_variant
        self.assertEqual(normalized_variant, [Variant('chr1', 10, 'G', '', self.genome)])

    def test_invalid_nucleotide(self):
        with self.assertRaisesRegex(ValueError, "Variant string was invalid:"):
            Variant.from_string('chr1:17:X:T', self.genome)

        with self.assertRaisesRegex(ValueError, "Variant string was invalid:"):
            Variant.from_string('chr1:16:A:X', self.genome)

    def test_compare(self):
        v = Variant("chr1", 2000, "A", "T", self.genome)
        self.assertEqual(v, Variant("chr1", 2000, "A", "T", self.genome))
        self.assertNotEqual(v, Variant("chr2", 2000, "A", "T", self.genome))
        self.assertNotEqual(v, Variant("chr1", 2001, "A", "T", self.genome))
        self.assertNotEqual(v, Variant("chr1", 2000, "C", "T", self.genome))
        self.assertNotEqual(v, Variant("chr1", 2000, "A", "C", self.genome))

    def test_sets(self):

        # If hash() or == for Variant is messed up, elements will be missing/repeated in sets.
        a = Variant("chr1", 2000, "A", "G", self.genome)
        b = Variant("chr1", 2000, "A", "G", self.genome)
        c = Variant("chr1", 2000, "C", "G", self.genome)
        d = Variant("chr1", 2000, "C", "G", self.genome)
        e = Variant("chr1", 2000, "A", "T", self.genome)
        f = Variant("chr1", 2000, "A", "T", self.genome)
        g = Variant("chr1", 2001, "A", "T", self.genome)
        h = Variant("chr1", 2001, "A", "T", self.genome)
        s1 = set([
            a,
            b,
            c,
            d,
            e,
            f,
            g,
            h,
        ])
        s2 = set([a, c, e, g])
        s3 = set([b, d, f, h])
        self.assertEqual(len(s1), 4)
        self.assertEqual(s1, s2)
        self.assertEqual(s1, s3)

    def test_spanning_error(self):
        interval = Interval('chr1', '+', 10, 12, self.genome, 10)
        with self.assertRaisesRegex(ValueError, "variants"):
            Variant.spanning(interval, interval)

        variants = [Variant('chr1', 10, 'N', 'A', self.genome)]
        with self.assertRaisesRegex(ValueError, "chromosome"):
            Variant.spanning(interval, Interval('chr2', interval.strand, interval.start, interval.end, interval.refg),
                             variants)
        with self.assertRaisesRegex(ValueError, "strand"):
            Variant.spanning(interval, interval.as_opposite_strand(), variants)
        with self.assertRaisesRegex(ValueError, "genome"):
            Variant.spanning(interval, Interval(interval.chrom, interval.strand, interval.start, interval.end, 'hg38.p12'),
                             variants)


    def test_spanning_fallback(self):
        x = Interval('chr1', '+', 10, 12, self.genome)
        y = x.shift(1)
        self.assertEqual(Variant.spanning(x, y), Interval.spanning(x, y))

    def test_spanning_point(self):
        variant = Variant('chr1', 9, 'N', 'NGTCGT', self.genome)
        interval = Interval('chr1', '+', 10, 12, self.genome, 10)

        start = interval.end5.with_anchor(None)
        middle = start.with_anchor(10, 1)
        end = start.with_anchor(10, 2)
        before = start.shift(-1)
        after = start.with_anchor(10, 3)
        far_after = start.shift(1)
        #  N  N  G  T  C  G  T  N  N  N
        # 8  9 10                11 12
        #       i  i  i
        #    b  s  m  e  a        f

        self.assertEqual(Variant.spanning(interval, start, [variant]), interval)
        self.assertEqual(Variant.spanning(interval, middle, [variant]), interval)
        self.assertEqual(Variant.spanning(interval, end, [variant]), interval)
        self.assertEqual(Variant.spanning(interval, before, [variant]), interval.expand(1, 0))
        self.assertEqual(Variant.spanning(interval, after, [variant]), interval.expand(0, 1))
        self.assertEqual(Variant.spanning(interval, far_after, [variant]), interval.expand(0, 4))

    def test_spanning_offsetted_point(self):
        variant = Variant('chr1', 9, 'N', 'NGTCGT', self.genome)
        interval = Interval('chr1', '-', 8, 10, self.genome, 10, 1)

        start = interval.end5.with_anchor(10, 1)
        middle = start.with_anchor(None)
        end = middle.shift(1)
        far_before = middle.shift(-1)
        before = start.with_anchor(10, 2)
        after = end.shift(1)
        #  N  N  G  T  C  G  T  N  N  N
        # 8  9 10                11 12
        #    i  i  i
        # a  e  m  s  b           f

        self.assertEqual(Variant.spanning(interval, start, [variant]), interval)
        self.assertEqual(Variant.spanning(interval, middle, [variant]), interval)
        self.assertEqual(Variant.spanning(interval, end, [variant]), interval)
        self.assertEqual(Variant.spanning(interval, far_before, [variant]), interval.expand(5, 0))
        self.assertEqual(Variant.spanning(interval, before, [variant]), interval.expand(1, 0))
        self.assertEqual(Variant.spanning(interval, after, [variant]), interval.expand(0, 1))

    def test_spanning_disjoint(self):
        variants = [Variant('chr1', 9, 'N', 'NGTCGT', self.genome), Variant('chr1', 10, 'N', 'NAAAAAAAA', self.genome)]
        left = Interval('chr1', '+', 10, 12, self.genome, 10, 1)
        right = Interval('chr1', '+', 11, 13, self.genome, 11, 1)
        #  N  N  G  T  C  G  T  N  A  A  A  A  A  A  A  A  N  N
        # 8  9 10                11                         12 13
        #          l  l  l           r  r  r
        self.assertEqual(Variant.spanning(left, right, variants), left.expand(0, 6))

        left = Interval('chr1', '+', 11, 13, self.genome, 11, 1)
        right = Interval('chr1', '+', 11, 13, self.genome, 11, 5)
        # N  A  A  A  A  A  A  A  A  N  N
        #  11                         12 13
        #      r  r  r     l  l  l
        self.assertEqual(Variant.spanning(left, right, variants), left.expand(0, 4))

    def test_spanning_overlapping(self):
        variants = [Variant('chr1', 9, 'N', 'NGTCGT', self.genome), Variant('chr1', 10, 'N', 'NAAAAA', self.genome)]
        left = Interval('chr1', '+', 10, 16, self.genome, 10, 1)
        right = Interval('chr1', '+', 11, 13, self.genome, 11, 1)
        #  N  N  G  T  C  G  T  N  A  A  A  A  A  N  N
        # 8  9 10                11                12 13
        #          l  l  l  l  l  l  l
        #                            r  r  r
        self.assertEqual(Variant.spanning(left, right, variants), left.expand(0, 2))

    def test_spanning_contains(self):
        variants = [Variant('chr1', 9, 'N', 'NGTCGT', self.genome), Variant('chr1', 10, 'N', 'NAAAAA', self.genome)]
        left = Interval('chr1', '+', 10, 21, self.genome, 10, 1)
        right = Interval('chr1', '+', 11, 13, self.genome, 11, 1)
        #  N  N  G  T  C  G  T  N  A  A  A  A  A  N  N
        # 8  9 10                11                12 13
        #          l  l  l  l  l  l  l  l  l  l  l  l
        #                            r  r  r
        self.assertEqual(Variant.spanning(left, right, variants), left)

    def test_spanning_deletion(self):
        variants = [Variant('chr1', 9, 'NNN', '', self.genome)]
        left = Interval('chr1', '+', 8, 8, self.genome)
        right = Interval('chr1', '+', 12, 13, self.genome, 13)
        #  N  -  -  -  N  N
        # 8  9 10 11 12 13 14
        #    r           r
        # l  l
        self.assertEqual(Variant.spanning(left, right, variants), right.expand(1, 0))

    def test_spanning_complex_substititon(self):
        variants = [Variant('chr1', 9, 'NNN', 'AT', self.genome)]
        left = Interval('chr1', '+', 9, 10, self.genome, 9)
        right = Interval('chr1', '+', 10, 11, self.genome)
        after = Interval('chr1', '+', 12, 13, self.genome)
        #  N  A  T  -  N  N  N
        # 8  9 10 11 12 13 14 15
        #    l  l
        #       r  r
        #             a  a
        self.assertEqual(Variant.spanning(left, right, variants), left.expand(0, 1))
        self.assertEqual(Variant.spanning(left, after, variants), left.expand(0, 2))

        unused_offset = Interval('chr1', '+', 9, 10, self.genome, 9, 1)
        self.assertEqual(Variant.spanning(unused_offset, after, variants), unused_offset.expand(0, 2))
