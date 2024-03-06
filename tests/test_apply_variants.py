# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from . import MiniGenomeDNA, MiniGenome
from genome_kit._apply_variants import apply_variants
from genome_kit._apply_variants import _apply_variants_left_anchor
from genome_kit._apply_variants import _apply_variants_right_anchor
from genome_kit._apply_variants import _apply_variants_no_anchor
from genome_kit import Interval
from genome_kit import Variant
import unittest

########################################################


class TestApplyVariantsRightAnchor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dna = MiniGenomeDNA('test_genome')

    @classmethod
    def tearDownClass(cls):
        del cls.dna

    def test_empty_variants(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        var_seq = _apply_variants_right_anchor(self.dna, [], interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_subst_internal(self):
        interval = Interval('chr1', '+', 10, 16, 'test_genome')
        variants = [Variant('chr1', 12, 'AC', 'GT', 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GTGTGT')

    def test_var_upstream(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        variants = [Variant('chr1', 3, "", "AAAAA", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_var_downstream(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        variants = [Variant('chr1', 30, "", "AAAAA", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_deletion_overlap_start(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 7, 'TACGT', "", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'CGACG')

    def test_indel_overlap_start(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 7, 'TACGT', "GG", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GGACG')

    def test_deletion_overlap_end(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 14, 'GTA', "", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'CGTAC')

    def test_deletion_interior(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 12, 'ACG', "", 'test_genome')]
        var_seq, alignment = _apply_variants_right_anchor(self.dna, variants, interval, True)
        self.assertEqual(var_seq, 'TACGT')

    def test_indel_overlap_end(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 13, 'CGTAC', "AAAA", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GTAAA')

    def test_double_variant(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 11, 'T', 'G', 'test_genome'), Variant('chr1', 13, 'CGTAC', "AAAA", 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GGAAA')

    def test_insertion_substitution_overlap_start(self):
        # Test a substitution followed by an insertion in the same variant

        # This simulates the effect of "chr10:76788941:GGAAAACCAG:GAAAAACCAAAA"
        # from Clinvar

        interval = Interval('chr1', '+', 15, 35, 'test_genome', 35)
        variants = [Variant('chr1', 11, 'TACGTACGT', 'GAAAAACCAAAA', 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq[:5], 'CAAAA')

        # This simulates the effect of "chr11:108141787:TTAGT:TGATACTA"
        # from Clinvar

        interval = Interval('chr1', '+', 16, 35, 'test_genome', 35)
        variants = [Variant('chr1', 16, 'ACGTA', 'TGATACTA', 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq[:5], 'TACTA')

    def test_deletion_substitution_overlap_start(self):
        # Test a substitution followed by a deletion

        interval = Interval('chr1', '+', 15, 35, 'test_genome', 35)
        variants = [Variant('chr1', 14, 'GTACGT', 'TTT', 'test_genome')]
        var_seq = _apply_variants_right_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq[:5], 'ACTTT')


class TestApplyVariantsLeftAnchor(unittest.TestCase):
    def setUp(self):
        self.genome = MiniGenome('test_genome')

    @classmethod
    def setUpClass(cls):
        cls.dna = MiniGenomeDNA('test_genome')

    @classmethod
    def tearDownClass(cls):
        del cls.dna

    def test_empty_variants(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        var_seq = _apply_variants_left_anchor(self.dna, [], interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_subst_internal(self):
        interval = Interval('chr1', '+', 10, 16, 'test_genome')
        variants = [Variant('chr1', 12, 'AC', 'GT', 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GTGTGT')

    def test_var_upstream(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        variants = [Variant('chr1', 3, "", "AAAAA", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_var_downstream(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        variants = [Variant('chr1', 30, "", "AAAAA", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_deletion_overlap_start(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 7, 'TACGT', "", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ACGTA')

    def test_indel_overlap_start(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 7, 'TACGT', "GG", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ACGTA')

    def test_deletion_overlap_end(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 14, 'GTA', "", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GTACC')

    def test_deletion_interior(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 12, 'ACG', "", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GTTAC')

    def test_indel_overlap_end(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 13, 'CGTAC', "AAAA", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GTAAA')

    def test_double_variant(self):
        interval = Interval('chr1', '+', 10, 15, 'test_genome')
        variants = [Variant('chr1', 11, 'T', 'G', 'test_genome'), Variant('chr1', 13, 'CGTAC', "AAAA", 'test_genome')]
        var_seq = _apply_variants_left_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, 'GGAAA')


class TestApplyVariantsNoAnchor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dna = MiniGenomeDNA('test_genome')

    @classmethod
    def tearDownClass(cls):
        del cls.dna

    def test_skipped_variant(self):
        interval = Interval('chr1', '+', 10, 20, 'test_genome')
        variants = [Variant('chr1', 5, 'A', '', 'test_genome')]
        var_seq = _apply_variants_no_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, self.dna(interval))

    def test_overlap_deletion(self):
        # Deletion overlaps the entire interval, must return empty sequence
        interval = Interval('chr1', '+', 10, 11, 'test_genome')
        variants = [Variant('chr1', 9, 'CGT', '', 'test_genome')]
        var_seq = _apply_variants_no_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, "")

        interval = Interval('chr1', '+', 10, 11, 'test_genome')
        variants = [Variant('chr1', 10, 'GTACG', '', 'test_genome')]
        var_seq = _apply_variants_no_anchor(self.dna, variants, interval)
        self.assertEqual(var_seq, "")


class TestApplyVariants(unittest.TestCase):
    def setUp(self):
        self.genome = MiniGenome('test_genome')

    @classmethod
    def setUpClass(cls):
        cls.dna = MiniGenomeDNA('test_genome')

    @classmethod
    def tearDownClass(cls):
        del cls.dna

    def test_empty_variants(self):
        # Forward strand
        interval = Interval('chr1', '+', 10, 20, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, [], interval)
        interval = Interval('chr1', '+', 10, 20, 'test_genome', 20)
        var_seq_left = apply_variants(self.dna, [], interval)
        interval = Interval('chr1', '+', 10, 20, 'test_genome', 15)
        var_seq_center = apply_variants(self.dna, [], interval)

        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, self.dna(interval))

        # Reverse strand
        interval = Interval('chr1', '-', 10, 20, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, [], interval)
        interval = Interval('chr1', '-', 10, 20, 'test_genome', 20)
        var_seq_left = apply_variants(self.dna, [], interval)
        interval = Interval('chr1', '-', 10, 20, 'test_genome', 15)
        var_seq_center = apply_variants(self.dna, [], interval)

        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, self.dna(interval))

    def test_subst_internal(self):
        variants = [Variant.from_string('chr1:13:AC:GT', self.genome)]
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 16)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, 'GTGTGT')

    def test_var_upstream(self):
        variants = [Variant.from_string("chr1:4:.:AAAAA", self.genome)]

        # Forward strand
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 16)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, self.dna(interval))

        # Reverse strand
        interval = Interval('chr1', '-', 10, 16, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '-', 10, 16, 'test_genome', 16)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '-', 10, 16, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, self.dna(interval))

    def test_var_downstream(self):
        variants = [Variant.from_string("chr1:31:-:AAAAA", self.genome)]

        # Forward strand
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 16)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 16, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, self.dna(interval))

        # Reverse strand
        interval = Interval('chr1', '-', 10, 16, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '-', 10, 16, 'test_genome', 16)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '-', 10, 16, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq_right, var_seq_left)
        self.assertEqual(var_seq_right, var_seq_center)
        self.assertEqual(var_seq_right, self.dna(interval))

    def test_deletion_overlap_start(self):
        variants = [Variant.from_string("chr1:8:TACGT:", self.genome)]

        # Forward strand
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 15)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq_right, 'ACGTA')
        self.assertEqual(var_seq_left, 'CGACG')
        self.assertEqual(var_seq_left, var_seq_center)

    def test_deletion_overlap_end(self):
        variants = [Variant.from_string("chr1:15:GTA:-", self.genome)]

        # Forward strand
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 15)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, 'GTACC')
        self.assertEqual(var_seq_left, 'CGTAC')
        self.assertEqual(var_seq_right, var_seq_center)

    def test_deletion_interior(self):
        variants = [Variant.from_string("chr1:13:ACG:.", self.genome)]

        interval = Interval('chr1', '+', 10, 15, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 15)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, 'GTTAC')
        self.assertEqual(var_seq_left, 'TACGT')
        self.assertEqual(var_seq_center, 'CGTTA')

    def test_indel_overlap_end(self):
        variants = [Variant.from_string("chr1:14:CGTAC:AAAA", self.genome)]

        # Forward strand
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 15)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, 'GTAAA')
        self.assertEqual(var_seq_left, 'GTAAA')
        self.assertEqual(var_seq_center, 'GTAAA')

    def test_double_variant(self):
        variants = [
            Variant.from_string("chr1:12:T:G", self.genome),
            Variant.from_string("chr1:14:CGTAC:AAAA", self.genome)
        ]

        interval = Interval('chr1', '+', 10, 15, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 15)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, 'GGAAA')
        self.assertEqual(var_seq_left, 'GGAAA')
        self.assertEqual(var_seq_center, 'GGAAA')

    def test_indel_overlap_entire_interval(self):
        variants = [Variant.from_string("chr1:9:ACGTACGT:", self.genome)]

        interval = Interval('chr1', '+', 10, 15, 'test_genome', 10)
        var_seq_right = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 15)
        var_seq_left = apply_variants(self.dna, variants, interval)
        interval = Interval('chr1', '+', 10, 15, 'test_genome', 13)
        var_seq_center = apply_variants(self.dna, variants, interval)

        self.assertEqual(var_seq_right, 'ACGTA')
        self.assertEqual(var_seq_left, 'NACGT')
        self.assertEqual(var_seq_center, 'CGTAC')

    def test_no_anchor(self):
        variants = [Variant.from_string('chr1:8:TACGTA:-', self.genome)]

        interval = Interval('chr1', '+', 20, 30, 'test_genome', None)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, self.dna(interval))

        interval = Interval('chr1', '+', 10, 15, 'test_genome', None)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'CG')

        interval = Interval('chr1', '+', 5, 10, 'test_genome', None)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'CG')

    def test_insertion_at_anchor(self):
        variant = [Variant.from_string('chr1:21::AGTT', self.genome)]

        interval = Interval('chr1', '+', 15, 25, 'test_genome', 20)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'TACGTAGTTA')

    def test_insertion_at_anchor_with_offset(self):
        variant = [Variant.from_string('chr1:21::AGTT', self.genome)]
        interval = Interval('chr1', '+', 15, 25, 'test_genome', 20, 2)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'CGTAGTTACG')

        variant = Variant.from_string("chr1:11::TTTTAGTTTT", self.genome)
        interval = Interval("chr1", "+", 10, 12, 'test_genome', 10, 4)
        sequence = apply_variants(self.dna, [variant], interval)
        self.assertEqual('AG', sequence)

    def test_anchor_outside_interval(self):
        # Deletion between interval end and anchor
        interval = Interval('chr1', '+', 10, 14, self.genome, 21)
        variant = Variant('chr1', 16, 'AC', '', self.genome)
        sequence = apply_variants(self.dna, [variant], interval)
        self.assertEqual(sequence, 'ACGT')

        # Deletion overlapping the end
        variant = [Variant.from_string('chr1:14:CG:', self.genome)]
        interval = Interval('chr1', '+', 10, 14, 'test_genome', 21)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'ACGT')

        # Insertion between interval end and anchor
        variant = [Variant.from_string('chr1:17::TTTTTTT', self.genome)]
        interval = Interval('chr1', '+', 10, 14, 'test_genome', 21)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'TTTT')

        # Deletion between the anchor and interval start
        variant = [Variant.from_string('chr1:7:GT:', self.genome)]
        interval = Interval('chr1', '+', 10, 14, 'test_genome', 5)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'ACGT')

        # Deletion overlapping the start
        variant = [Variant.from_string('chr1:10:CG:', self.genome)]
        interval = Interval('chr1', '+', 10, 14, 'test_genome', 5)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'ACGT')

        # Insertion between anchor and interval start
        variant = [Variant.from_string('chr1:9::TTTTTTT', self.genome)]
        interval = Interval('chr1', '+', 10, 14, 'test_genome', 5)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual(sequence, 'TTTT')

    def test_insertion_at_interval_ends(self):
        # Test behaviour when insertions occur at the start or
        # end of an interval, with or without an anchor

        # Insertion at start
        variant = [Variant.from_string('chr1:10::TTT', self.genome)]

        # Insertion will be inside interval
        # (implies anchor_offset=0)
        interval = Interval('chr1', '+', 9, 19, 'test_genome', 9)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual('TTTCGTACGT', sequence)

        # Insertion will be outside interval
        interval = Interval('chr1', '+', 9, 19, 'test_genome', 9, 3)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual('CGTACGTACG', sequence)

        # Insertion at end
        variant = [Variant.from_string('chr1:20::TTT', self.genome)]

        # Insertion will be outside interval
        # (implies anchor_offset=0)
        interval = Interval('chr1', '+', 9, 19, 'test_genome', 19)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual('CGTACGTACG', sequence)

        # Insertion will be inside interval
        interval = Interval('chr1', '+', 9, 19, 'test_genome', 19, 3)
        sequence = apply_variants(self.dna, variant, interval)
        self.assertEqual('ACGTACGTTT', sequence)

    def test_reference_alignment_no_anchor(self):
        genome37 = MiniGenome('test_genome')
        variants = [Variant.from_string("chr1:11:G:-", genome37)]
        interval = Interval('chr1', '+', 5, 15, genome37)

        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]

        self.assertEqual([0, 1, 2, 3, 4, 6, 7, 8, 9], reference_alignment)

        variants = [Variant.from_string("chr1:11:GTA:TT", genome37)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 8, 9], reference_alignment)

        variants = [Variant.from_string("chr1:11::TT", genome37)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual([0, 1, 2, 3, 4, (5, 0), (5, 1), 5, 6, 7, 8, 9], reference_alignment)

        variants = [Variant.from_string("chr1:11:G:TT", genome37)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual([0, 1, 2, 3, 4, 5, (6, 0), 6, 7, 8, 9], reference_alignment)

    def test_reference_alignment_interior_anchor(self):
        genome37 = MiniGenome('test_genome')
        interval = Interval('chr1', '+', 5, 15, genome37, 10)

        # Test insertion
        variants = [Variant.from_string("chr1:11::TT", self.genome)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]

        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, (5, 0), (5, 1), 5, 6, 7])

        # Test substitution
        variants = [Variant.from_string("chr1:11:GTA:ATG", genome37)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]

        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Test insertion with anchor_offset
        interval = Interval('chr1', '+', 5, 15, genome37, 10, 2)
        variants = [Variant.from_string("chr1:11::TTATT", self.genome)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]

        self.assertEqual(reference_alignment, [2, 3, 4, (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), 5, 6])

        # Test deletion
        interval = Interval('chr1', '+', 5, 15, genome37, 10)
        variants = [Variant.from_string("chr1:11:GTA:", genome37)]
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, 8, 9, 10, 11, 12])

    def test_reference_alignment_other_anchor_cases(self):
        genome37 = MiniGenome('test_genome')
        variants = [Variant.from_string("chr1:11:GTA:ATG", genome37)]

        # Test anchor==start
        interval = Interval('chr1', '+', 5, 15, genome37, 5)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]

        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Test anchor==start
        interval = Interval('chr1', '+', 5, 15, genome37, 15)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Test anchor > end
        interval = Interval('chr1', '+', 5, 15, genome37, 20)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Test anchor < start with substitution
        interval = Interval('chr1', '+', 5, 15, genome37, 0)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

        # Test anchor < start with substitution
        variants = [Variant.from_string("chr1:11::TT", self.genome)]
        interval = Interval('chr1', '+', 5, 15, genome37, 0)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [0, 1, 2, 3, 4, (5, 0), (5, 1), 5, 6, 7])

        # Test anchor right of deletion
        variants = [Variant.from_string("chr1:11:GT:", self.genome)]
        interval = Interval('chr1', '+', 5, 15, genome37, 11)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [-1, 0, 1, 2, 3, 4, 6, 7, 8, 9])

        # Test anchor right of insertion
        variants = [Variant.from_string("chr1:11::TT", self.genome)]
        interval = Interval('chr1', '+', 5, 15, genome37, 15)
        reference_alignment = apply_variants(genome37.dna, variants, interval, reference_alignment=True)[1]
        self.assertEqual(reference_alignment, [2, 3, 4, (5, 0), (5, 1), 5, 6, 7, 8, 9])

    def test_reference_alignment_exception(self):
        genome37 = MiniGenome('test_genome')
        variants = [Variant.from_string("chr1:11:GTA:ATG", genome37)]
        interval = Interval('chr1', '-', 5, 15, genome37, 5)

        with self.assertRaises(ValueError):
            apply_variants(genome37.dna, variants, interval, reference_alignment=True)

    def test_variant_on_other_chromosome(self):
        """Tests that variants are only applied when they are on the same
        chromosome as the interval.
        """

        genome37 = MiniGenome('test_genome')
        variants = [Variant.from_string("chr1:11:GTA:ATG", genome37)]
        interval = Interval('chr2', '+', 10, 20, genome37)
        wt_dna = genome37.dna(interval)

        mt_dna = apply_variants(genome37.dna, variants, interval)
        self.assertEqual(wt_dna, mt_dna)

    def test_insertion_and_substitution(self):
        variants = [Variant.from_string('chr1:14::T', self.genome), Variant.from_string('chr1:15:G:A', self.genome)]
        self.assertEqual(self.genome.dna(Interval('chr1', '+', 11, 17, 'test_genome')), 'TACGTA')

        # no anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome')
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ATCAT')

        # right anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome', 16)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'TCAT')

        # at anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome', 13)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ATCA')

    def test_insertion_and_deletion(self):
        variants = [Variant.from_string('chr1:14::T', self.genome), Variant.from_string('chr1:15:G:', self.genome)]
        self.assertEqual(self.genome.dna(Interval('chr1', '+', 11, 17, 'test_genome')), 'TACGTA')

        # no anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome')
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ATCT')

        # right anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome', 16)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ATCT')

        # at anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome', 13)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'ATCT')

    def test_deletion_and_substitution(self):
        variants = [Variant.from_string('chr1:14:C:', self.genome), Variant.from_string('chr1:15:G:A', self.genome)]
        self.assertEqual(self.genome.dna(Interval('chr1', '+', 11, 17, 'test_genome')), 'TACGTA')

        # no anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome')
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'AAT')

        # right anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome', 16)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'TAAT')

        # at anchor
        interval = Interval('chr1', '+', 12, 16, 'test_genome', 13)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'AATA')

    def test_substitution_in_extended_interval_from_deletion(self):
        self.assertEqual(self.genome.dna(Interval('chr1', '+', 11, 17, 'test_genome')), 'TACGTA')

        # right anchor
        variants = [Variant.from_string('chr1:14:C:A', self.genome), Variant.from_string('chr1:15:G:', self.genome)]
        interval = Interval('chr1', '+', 14, 15, 'test_genome', 15)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'A')

        # left anchor
        variants = [Variant.from_string('chr1:14:C:', self.genome), Variant.from_string('chr1:15:G:A', self.genome)]
        interval = Interval('chr1', '+', 14, 15, 'test_genome', 14)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'A')

    def test_deletion_and_insertion(self):
        variants = [Variant.from_string('chr1:14:C:', self.genome), Variant.from_string('chr1:15::T', self.genome)]
        self.assertEqual(self.genome.dna(Interval('chr1', '+', 11, 17, 'test_genome')), 'TACGTA')

        # no anchor
        interval = Interval('chr1', '+', 13, 16, 'test_genome')
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'TGT')

        # right anchor
        interval = Interval('chr1', '+', 13, 16, 'test_genome', 16)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'TGT')

        # at anchor
        interval = Interval('chr1', '+', 13, 16, 'test_genome', 13)
        var_seq = apply_variants(self.dna, variants, interval)
        self.assertEqual(var_seq, 'TGT')
