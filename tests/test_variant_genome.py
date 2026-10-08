# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import gzip
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch
from genome_kit import DataManager
from genome_kit import VariantGenome
from genome_kit import Interval
from genome_kit import Variant
from genome_kit import gk_data
from genome_kit.data_manager import GKDataFileNotFoundError
from . import MiniGenome, get_test_file


class LocalDataManager(DataManager):
    """Stores uploads in a local directory, serving other files from the test data."""

    def get_file(self, filename):
        path = os.path.join(self.data_dir, filename)
        if os.path.exists(path):
            return path
        if filename.startswith(VariantGenome._FILENAME_FORMAT.split("{")[0]):
            raise GKDataFileNotFoundError(filename, "not found")
        return get_test_file(filename)

    def upload_file(self, filepath, filename, metadata=None):
        shutil.copyfile(filepath, os.path.join(self.data_dir, filename))


class TestVariantGenome(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        cls.genome = MiniGenome('test_genome')

    @classmethod
    def tearDownClass(cls):
        del cls.genome

    def test_as_genome(self):
        # Test that Genome attributes and methods accessible through VariantGenome
        # even if they were not explicitly overridden on VariantGenome. Use a Genome
        # instance that has annotations (unlike the current self.genome)
        anno_genome = MiniGenome("gencode.v29")
        variant = VariantGenome(anno_genome, Variant.from_string('chr2:20:A:T', anno_genome))
        self.assertEqual(variant.reference_genome, anno_genome.reference_genome)
        self.assertEqual(variant.data_dir, anno_genome.data_dir)
        self.assertIs(variant.annotation, anno_genome.annotation)
        self.assertIs(variant.genes, anno_genome.genes)
        self.assertIs(variant.exons, anno_genome.exons)

    def test_invalid_variant(self):
        with self.assertRaises(ValueError):
            VariantGenome(self.genome, Variant.from_string('chr99:1:A:G', self.genome))

        with self.assertRaises(ValueError):
            VariantGenome(self.genome, Variant.from_string('chr99-1-A-G', self.genome))

        with self.assertRaises(ValueError):
            VariantGenome(self.genome, Variant.from_string('chr1:foo:A:G', self.genome))

        with self.assertRaises(ValueError):
            VariantGenome(self.genome, Variant.from_string('chr1:1:A:G', self.genome))

    def test_unicode_variant(self):
        variant = VariantGenome(self.genome, Variant.from_string(u'chr1:15:G:T', self.genome))
        self.assertIsInstance(variant.__repr__(), str)
        self.assertIsInstance(variant.__str__(), str)
        for genome in (variant, self.genome):
            sequence = variant.dna(Interval('chr1', '+', 10, 20, genome, None))
            self.assertEqual(sequence[5], 'T')

        interval = Interval('chr1', '-', 10, 20, variant, None)
        sequence = variant.dna(interval)
        self.assertEqual(sequence[5], 'A')

    def test_single_snv(self):
        variant = VariantGenome(self.genome, Variant.from_string('chr1:15:G:T', self.genome))
        self.assertIsInstance(variant.__repr__(), str)
        self.assertIsInstance(variant.__str__(), str)
        for genome in (variant, self.genome):
            sequence = variant.dna(Interval('chr1', '+', 10, 20, genome, None))
            self.assertEqual(sequence[5], 'T')

        interval = Interval('chr1', '-', 10, 20, variant, None)
        sequence = variant.dna(interval)
        self.assertEqual(sequence[5], 'A')

    def test_deletion(self):
        variant = VariantGenome(self.genome, Variant.from_string('chr1:15:G:-', self.genome))
        interval = Interval('chr1', '+', 10, 20, variant, None)
        sequence = variant.dna(interval)
        self.assertEqual(len(sequence), len(interval) - 1)
        self.assertEqual(sequence, 'GTACTACGT')

    def test_clinvar_deletion_format(self):
        variant_clinvar = VariantGenome(self.genome, Variant.from_string('chr2:11:GT:G', self.genome))
        variant_non_clinvar = VariantGenome(self.genome, Variant.from_string('chr2:12:T:', self.genome))

        interval = Interval('chr1', '+', 10, 20, variant_clinvar, 'start')

        self.assertEqual(variant_clinvar.dna(interval), variant_non_clinvar.dna(interval))

    def test_find_motif(self):
        # For reference, test_genome's chr2 is 40bp long and
        # contains the following sequence:
        #
        #          0         1         2         3         4
        #          0123456789012345678901234567890123456789
        #   chr2 = AACCTTTTACGTAAACCCGGGTTTACCGGAAATGGATTAA
        #
        variant = Variant.from_string("chr2:11:G:A", self.genome)
        variant_genome = VariantGenome(self.genome, variant)

        # Test SNV on the forward strand
        interval = Interval('chr1', '+', 0, 40, self.genome)
        motif = 'GT'
        motif_hits = variant_genome.find_motif(interval, motif)

        self.assertEqual(motif, variant_genome.dna(motif_hits[0].expand(0, 2)))

        # Test a deletion on the forward strand
        variant = Variant.from_string("chr2:11:G:-", self.genome)
        variant_genome = VariantGenome(self.genome, variant)
        interval = Interval('chr2', '+', 5, 15, self.genome)

        # Deletion doesn't change coordinate of this motif
        motif = 'TAC'
        motif_hits = variant_genome.find_motif(interval, motif)
        dna = self.genome.dna(motif_hits[0].expand(0, 3))
        self.assertEqual(motif, dna)

        # Deletion will change the position of this motif
        motif = 'TAA'
        motif_hit_wt = self.genome.find_motif(interval, motif)
        motif_hit_mt = variant_genome.find_motif(interval, motif)

        self.assertEqual(motif_hit_wt[0].start, motif_hit_mt[0].start)
        self.assertEqual(motif_hit_wt[0].end, motif_hit_mt[0].end)

        # Motif is contained in insertion
        variant = Variant.from_string("chr2:11::TTTTAGTTTT", self.genome)
        variant_genome = VariantGenome(self.genome, variant)
        motif = 'AG'
        motif_hits = variant_genome.find_motif(interval, motif)
        motif_interval = motif_hits[0].expand(0, 2)
        dna = variant_genome.dna(motif_interval)
        self.assertEqual(motif, dna)

        # Try the same on reverse strand
        interval = Interval('chr2', '-', 5, 15, self.genome)
        motif = 'CT'
        motif_hits = variant_genome.find_motif(interval, motif)
        motif_interval = motif_hits[0].expand(0, 2)
        dna = variant_genome.dna(motif_interval)
        self.assertEqual(motif, dna)

        # Test string arguments
        interval = Interval('chr2', '+', 0, 39, self.genome)
        motif = 'AG'

        motif_hits = variant_genome.find_motif(interval, motif, '5p')
        self.assertEqual(motif, variant_genome.dna(motif_hits[0].expand(0, 2)))

        motif_hits = variant_genome.find_motif(interval, motif, '3p')
        self.assertEqual(motif, variant_genome.dna(motif_hits[0].expand(2, 0)))

        # Test match_position in range
        with self.assertRaises(ValueError):
            variant_genome.find_motif(interval, motif, 3)

        # Test motif occurs at end of interval with match_position=len(motif)
        variant = Variant.from_string("chr2:1:A:T", self.genome)
        variant_genome = VariantGenome(self.genome, variant)
        interval = Interval('chr2', '+', 0, 12, self.genome)
        motif = 'CGT'
        motif_hits = variant_genome.find_motif(interval, motif, match_position=len(motif))
        dna = variant_genome.dna(motif_hits[0].expand(len(motif), 0))
        self.assertEqual(motif, dna)

        # Test overlapping and non-overlapping matching
        motif = 'TT'
        variant = Variant.from_string('chr2:20:G:A', self.genome)  # Use variant outside interval
        variant_genome = VariantGenome(self.genome, variant)
        interval = Interval('chr2', '+', 0, 10, self.genome)
        motif_hits = variant_genome.find_motif(interval, motif, find_overlapping_motifs=False)
        self.assertEqual([Interval("chr2", "+", 4, 4, self.genome, 4), Interval("chr2", "+", 6, 6, self.genome, 6)], motif_hits)

        motif_hits = variant_genome.find_motif(interval, motif, find_overlapping_motifs=True)
        self.assertEqual([
            Interval("chr2", "+", 4, 4, self.genome, 4),
            Interval("chr2", "+", 5, 5, self.genome, 5),
            Interval("chr2", "+", 6, 6, self.genome, 6)
        ], motif_hits)

        # Test 3p at an insertion
        motif = 'TA'
        variant = Variant.from_string('chr2:9::A', self.genome)
        variant_genome = VariantGenome(self.genome, variant)
        interval = Interval('chr2', '+', 0, 10, self.genome)
        motif_hits = variant_genome.find_motif(interval, motif, '3p')
        dna = variant_genome.dna(motif_hits[0].expand(len(motif), 0))
        self.assertEqual(motif, dna)

    def test_reference_variant_genome_match(self):
        # This runs the tests from TestGenome.test_find_motif in both reference and variant
        # genome and ensures that both return the same results

        genome37 = MiniGenome("hg19")
        variant_genome = VariantGenome(genome37, [])
        motif = 'AG'

        # Test on forward strand without position_offset
        interval_forward = Interval('chr2', '+', 100, 200, genome37)
        motif_hits_ref = genome37.find_motif(interval_forward, motif)
        motif_hits_var = variant_genome.find_motif(interval_forward, motif)

        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        # Test on forward strand with position offset
        motif_hits_ref = genome37.find_motif(interval_forward, motif, len(motif))
        motif_hits_var = variant_genome.find_motif(interval_forward, motif, len(motif))

        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        # Test on reverse strand without position offset
        interval_reverse = Interval('chr2', '-', 100, 200, genome37)

        motif_hits_ref = genome37.find_motif(interval_reverse, motif)
        motif_hits_var = variant_genome.find_motif(interval_reverse, motif)

        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        # Test on reverse strand with position offset
        motif_hits_ref = genome37.find_motif(interval_reverse, motif, len(motif))
        motif_hits_var = variant_genome.find_motif(interval_reverse, motif, len(motif))

        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        # Test string arguments
        interval = Interval('chr2', '+', 0, 100, genome37)
        motif = 'AA'

        motif_hits_ref = genome37.find_motif(interval, motif, '5p')
        motif_hits_var = variant_genome.find_motif(interval, motif, '5p')
        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        motif_hits_ref = genome37.find_motif(interval, motif, '3p')
        motif_hits_var = variant_genome.find_motif(interval, motif, '3p')
        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        # Test match_position in range
        with self.assertRaises(ValueError):
            genome37.find_motif(interval, motif, 3)
            variant_genome.find_motif(interval, motif, 3)

        # Test overlapping and non-overlapping matching
        motif = 'TT'
        interval = Interval('chr2', '+', 0, 10, genome37)
        motif_hits_ref = genome37.find_motif(interval, motif, find_overlapping_motifs=False)
        motif_hits_var = variant_genome.find_motif(interval, motif, find_overlapping_motifs=False)
        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

        motif_hits_ref = genome37.find_motif(interval, motif, find_overlapping_motifs=True)
        motif_hits_var = variant_genome.find_motif(interval, motif, find_overlapping_motifs=True)
        self.assertEqual(motif_hits_ref, [x.with_anchor(None) for x in motif_hits_var])

    def test_non_matching_genome(self):
        genome38 = MiniGenome('hg38.p12')

        with self.assertRaisesRegex(ValueError, "genome doesn't match the reference genome."):
            VariantGenome(self.genome, Variant.from_string("chr2:11:G:A", genome38))

    def test_string_variant(self):
        """Test that trying to initialize with a variant string raises an
        error"""

        with self.assertRaises(TypeError):
            VariantGenome(self.genome, "chr1:15:G:T")

        with self.assertRaises(TypeError):
            VariantGenome(self.genome, ["chr1:15:G:T"])

    def _variant(self, chrom, start, alt, ref_len=1, genome=None):
        genome = genome or self.genome
        ref = genome.dna(Interval(chrom, "+", start, start + ref_len, genome))
        return Variant(chrom, start, ref, alt, genome)

    def test_eq(self):
        # Same start on different chromosomes, in different orders
        a = self._variant("chr1", 10, "T")
        b = self._variant("chr2", 10, "C")
        self.assertEqual(VariantGenome(self.genome, [a, b]), VariantGenome(self.genome, [b, a]))

        # Un-normalized variant equals its normalized form
        ref = self.genome.dna(Interval("chr1", "+", 10, 12, self.genome))
        # G-padded = "11:GT:GC", non-padded = "12:T:C"
        padded = Variant("chr1", 10, ref, ref[0] + "C", self.genome)
        non_padded = self._variant("chr1", 11, "C")
        self.assertEqual(VariantGenome(self.genome, padded), VariantGenome(self.genome, non_padded))

        # different variants
        self.assertNotEqual(VariantGenome(self.genome, a), VariantGenome(self.genome, b))
        # different genome
        self.assertNotEqual(VariantGenome(self.genome, a), VariantGenome(MiniGenome("hg19"), self._variant("chr1", 10, "T", genome=MiniGenome("hg19"))))

    def test_id(self):
        a = self._variant("chr1", 10, "T")
        b = self._variant("chr2", 10, "C")
        vg = VariantGenome(self.genome, [a, b])
        self.assertRegex(vg.id, "^[0-9a-f]{64}$")
        self.assertEqual(vg.id, VariantGenome(self.genome, [b, a]).id)
        self.assertNotEqual(vg.id, VariantGenome(self.genome, a).id)
        self.assertNotEqual(vg.id, VariantGenome(self.genome, []).id)

    def test_save_load(self):
        variants = [self._variant("chr1", 10, "T"), self._variant("chr2", 10, "C"), self._variant("chr1", 20, "", ref_len=2)]
        vg = VariantGenome(self.genome, variants)
        with tempfile.TemporaryDirectory() as tmpdir, patch.object(gk_data, "data_manager", LocalDataManager(tmpdir)):
            vg_id = vg.save()
            self.assertEqual(vg_id, vg.id)
            with open(os.path.join(tmpdir, VariantGenome._FILENAME_FORMAT.format(vg_id)), "rb") as f:
                saved = f.read()
            with patch("time.time", return_value=0.0):
                VariantGenome(self.genome, variants[::-1]).save()
            with open(os.path.join(tmpdir, VariantGenome._FILENAME_FORMAT.format(vg_id)), "rb") as f:
                self.assertEqual(f.read(), saved)
            loaded = VariantGenome.load(vg_id)
            self.assertEqual(loaded, vg)
            self.assertEqual(loaded.id, vg_id)
            self.assertEqual(loaded.dna(Interval("chr1", "+", 5, 25, loaded)), vg.dna(Interval("chr1", "+", 5, 25, vg)))

            with self.assertRaises(GKDataFileNotFoundError):
                VariantGenome.load("0" * 64)

            for bad_id in ["", vg_id[:-1], vg_id.upper(), "../" + vg_id, vg_id + "\n"]:
                with self.assertRaisesRegex(ValueError, "Invalid variant genome id"):
                    VariantGenome.load(bad_id)

            # Stored content that does not match its id is rejected
            with gzip.open(os.path.join(tmpdir, VariantGenome._FILENAME_FORMAT.format(vg_id)), "wt") as f:
                f.write(self.genome.config + "\n")
            with self.assertRaisesRegex(ValueError, "does not match id"):
                VariantGenome.load(vg_id)

    def test_allow_outside_chromsome(self):
        variant = VariantGenome(
            self.genome, Variant.from_string("chr1:15:G:T", self.genome)
        )
        sequence = variant.dna(variant.interval("chr1", "+", -10, 20))
        self.assertEqual(len(sequence), 30)
        self.assertEqual(sequence[0], "N")
        self.assertEqual(sequence[-5], "T")
