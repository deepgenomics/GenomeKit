# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import gc
import unittest

from genome_kit import Genome, Interval, Variant

from . import MiniGenome

TEST_REFGS = ["hg19", "hg19.p13.plusMT", "hg38.p12"]

########################################################


class TestGenome(unittest.TestCase):
    def test_instantiate(self):
        # Create and re-create several Genome instances.
        # Make sure they're deallocated properly.
        for i in range(3):
            g = MiniGenome()  # Default annotation, whatever that is
            self.assertIsInstance(g.data_dir, str)
            self.assertIsInstance(g.refg, str)
            self.assertTrue("hg19" in g.refg)
            self.assertIsInstance(g.__repr__(), str)
            self.assertEqual(g._chromosome_sizes, None)

            del g
            gc.collect()
            self.assertEqual(gc.garbage, [])

        MiniGenome("gencode.v29")
        MiniGenome("gencode.v29lift37")

    def test_attributes(self):
        # Load a specific annotation and check attributes
        g = MiniGenome()
        self.assertIsInstance(g.data_dir, str)
        self.assertIsInstance(g.config, str)
        self.assertIsInstance(g.reference_genome, str)
        self.assertIsInstance(g.__repr__(), str)
        self.assertEqual(g._chromosome_sizes, None)

        with self.assertRaises(TypeError):
            g.data_dir = "gencode.v20lift30"
        with self.assertRaises(TypeError):
            g.reference_genome = "h15"
        with self.assertRaises(Exception):
            g.dna = "test"
        with self.assertRaises(Exception):
            g.annotation = "test"

    def test_equal(self):
        a = MiniGenome("hg19")
        b = MiniGenome("hg19")
        c = MiniGenome("hg38.p12")
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)
        self.assertNotEqual(b, c)

    def test_incompatible_reference_genome(self):
        genome37 = MiniGenome("hg19")
        genome38 = MiniGenome("hg38.p12")
        # Genome to avoid creating mini metadata to pass the interval creation

        # Test with string reference_genome identifiers
        i19 = Interval("chr2", "+", 100, 200, genome37)
        i38 = Interval("chr2", "+", 100, 200, genome38)
        with self.assertRaises(ValueError):
            genome37.dna(i38)
        with self.assertRaises(ValueError):
            genome38.dna(i19)
        self.assertEqual(genome37.reference_genome, i19.reference_genome)
        self.assertEqual(genome38.reference_genome, i38.reference_genome)

        # Test with actual Genome instance as reference_genome
        i19 = Interval("chr2", "+", 100, 200, genome37)
        i38 = Interval("chr2", "+", 100, 200, genome38)
        with self.assertRaises(ValueError):
            genome37.dna(i38)
        with self.assertRaises(ValueError):
            genome38.dna(i19)
        self.assertEqual(genome37.reference_genome, i19.reference_genome)
        self.assertEqual(genome38.reference_genome, i38.reference_genome)

    def test_find_motif(self):
        genome37 = MiniGenome("test_genome")
        motif = 'AG'

        # Test on forward strand without position_offset
        interval_forward = Interval('chr2', '+', 5, 30, genome37)
        motif_hits = genome37.find_motif(interval_forward, motif)

        for hit in motif_hits:
            dna = genome37.dna(hit.expand(0, len(motif)))
            self.assertEqual(motif, dna)

        # Test on forward strand with position offset
        motif_hits = genome37.find_motif(interval_forward, motif, len(motif))

        for hit in motif_hits:
            dna = genome37.dna(hit.expand(len(motif), 0))
            self.assertEqual(motif, dna)

        # Test on reverse strand without position offset
        interval_reverse = Interval('chr2', '-', 5, 30, genome37)

        motif_hits = genome37.find_motif(interval_reverse, motif)

        for hit in motif_hits:
            dna = genome37.dna(hit.expand(0, len(motif)))
            self.assertEqual(motif, dna)

        # Test on reverse strand with position offset
        motif_hits = genome37.find_motif(interval_reverse, motif, len(motif))

        for hit in motif_hits:
            dna = genome37.dna(hit.expand(len(motif), 0))
            self.assertEqual(motif, dna)

        # Test string arguments
        interval = Interval('chr2', '+', 5, 35, genome37)
        motif = 'AA'

        motif_hits = genome37.find_motif(interval, motif, '5p')
        self.assertEqual(motif, genome37.dna(motif_hits[0].expand(0, 2)))

        motif_hits = genome37.find_motif(interval, motif, '3p')
        self.assertEqual(motif, genome37.dna(motif_hits[0].expand(2, 0)))

        # Test match_position in range
        with self.assertRaises(ValueError):
            genome37.find_motif(interval, motif, 3)

        # Test overlapping and non-overlapping matching
        motif = 'TT'
        interval = Interval('chr2', '+', 0, 10, genome37)
        motif_hits = genome37.find_motif(interval, motif, find_overlapping_motifs=False)
        self.assertEqual([Interval("chr2", "+", 4, 4, genome37), Interval("chr2", "+", 6, 6, genome37)], motif_hits)

        motif_hits = genome37.find_motif(interval, motif, find_overlapping_motifs=True)
        self.assertEqual([
            Interval("chr2", "+", 4, 4, genome37),
            Interval("chr2", "+", 5, 5, genome37),
            Interval("chr2", "+", 6, 6, genome37)
        ], motif_hits)

    def test_chromosome_size(self):
        for refg in TEST_REFGS:
            genome = Genome(refg)

            # when known chromosome
            chromosome = "chr1"
            self.assertIsInstance(genome.chromosome_size(chromosome), int)

            # when unknown chromosome
            chromosome = "chrXKASD"

            with self.assertRaises(ValueError):
                genome.chromosome_size(chromosome)

    def test__chrom_sizes(self):
        for refg in TEST_REFGS:
            # when self._chromosome_sizes is not empty
            genome = Genome(refg)
            genome._chromosome_sizes = {"not": "empty"}

            # it returns the already loaded chromosome sizes
            self.assertEqual(genome._chrom_sizes(), genome._chromosome_sizes)

            # when self._chromosome_sizes is empty
            genome._chromosome_sizes = {}

            # it deserializes the chromosome sizes from the disk into a str/int dict
            self.assertGreater(len(genome._chrom_sizes().keys()), 0)
            for key in genome._chrom_sizes().keys():
                self.assertIsInstance(key, str)

            for value in genome._chrom_sizes().values():
                self.assertIsInstance(value, int)

    def test_interval(self):
        for refg in TEST_REFGS:
            genome = Genome(refg)

            chromosome, strand, start, end = "chr1", "+", 10000, 10010
            interval = genome.interval(chromosome, strand, start, end)

            # it returns an interval with the same refg as the genome object it's being called on
            self.assertEqual(interval, Interval(chromosome, strand, start, end, genome.refg))

    def test_variant(self):
        genome = MiniGenome("hg19")

        variant_string = "chr2:1,10:A:G"
        variant = genome.variant(variant_string)

        # it returns a variant with the same refg as the genome object it's being called on
        self.assertIsInstance(variant, Variant)
        self.assertEqual(variant.reference_genome, genome.reference_genome)
        self.assertEqual(variant, Variant.from_string(variant_string, genome))

    def test_variant_dna(self):
        genome = MiniGenome("hg19")

        variant_string = "chr2:1,10:A:G"
        variant = genome.variant(variant_string)
        interval = Interval("chr2", "+", 109, 110, genome)

        self.assertEqual(genome.variant_dna(interval, variant), "G")

        interval = Interval("chr2", "+", -10, 110, genome)

        self.assertEqual(genome.variant_dna(interval, variant)[-1], "G")

    def test_appris_simple(self):
        genome = MiniGenome("gencode.v29")
        self.assertTrue(set(genome.appris_transcripts()).issubset(set(genome.transcripts)))

    def test_appris_transcripts_unavailable(self):
        genome = MiniGenome("gencode.v29lift37")
        with self.assertRaisesRegex(RuntimeError, "APPRIS not available"):
            genome.appris_transcripts()

    def test_appris_principality_exception(self):
        # Test for exception raising
        genome = MiniGenome("hg19")
        with self.assertRaisesRegex(ValueError, "Genome is not a registered annotation"):
            genome.appris_principality('ENST00000233331')

    def test_appris_transcripts_shallow_copy(self):
        """GK uses a list of transcript indices for quick queries of TranscriptTable, but should be
         careful not to return this list or mutate it somehow. The expected behaviour is to return a shallow
         copy of the list"""
        genome = MiniGenome("gencode.v29")

        transcripts = genome.appris_transcripts('ENSG00000274049')
        self.assertEqual(transcripts, genome.appris_transcripts('ENSG00000274049'))

        del transcripts[0]
        self.assertEqual(len(genome.appris_transcripts('ENSG00000274049')), len(transcripts) + 1)

        transcripts = genome.appris_transcripts('ENSG00000274049')
        transcripts[0] = None
        self.assertNotEqual(transcripts, genome.appris_transcripts('ENSG00000274049'))

    def test_appris_principality_gencode29(self):
        # Manually curated set of APPRIS principalities based on APPRIS gencode29 (ensembl 94)
        answers = [('ENSG00000115274', 'ENST00000233331', 0),
                   ('ENSG00000115274', 'ENST00000431187', None),
                   ('ENSG00000115274', 'ENST00000494986', None),
                   ('ENSG00000115274', 'ENST00000409917', None),
                   ('ENSG00000274049', 'ENST00000452361', 4),
                   ('ENSG00000274049', 'ENST00000441673', 6),
                   ('ENSG00000239779', 'ENST00000233615', 1),
                   ('ENSG00000239779', 'ENST00000393972', 5),
                   ('ENSG00000239779', 'ENST00000474185', None),
                   ('ENSG00000239779', 'ENST00000466835', None),
                   ('ENSG00000239779', 'ENST00000409737', 5),
                   ('ENSG00000239779', 'ENST00000428943', None),
                   ('ENSG00000239779', 'ENST00000494741', None),
                   ('ENSG00000239779', 'ENST00000473467', None),
                   ('ENSG00000115275', 'ENST00000233616', 2),
                   ('ENSG00000115275', 'ENST00000452063', 6),
                   ('ENSG00000115275', 'ENST00000409065', None),
                   ('ENSG00000115275', 'ENST00000462189', None),
                   ('ENSG00000115275', 'ENST00000448666', 6)]  # yapf: disable

        # Test on a small subset of gencode.v29 annotations
        genome = MiniGenome("gencode.v29")
        for g, t, a in answers:
            # Test for correct APPRIS scores using strings
            self.assertEqual(genome.appris_principality(t), a, t)
            # Test appris_transcripts()
            self.assertEqual(genome.appris_principality(genome.transcripts[t]), a)

        # Test appris_transcripts
        transcripts = genome.appris_transcripts('ENSG00000274049')
        answer = [("ENST00000452361.5", 4), ("ENST00000441673.2", 6)]
        for t, a in zip(transcripts, answer):
            self.assertEqual((t.id, genome.appris_principality(t)), a)

        transcripts = genome.appris_transcripts('ENSG00000115275')
        answer = [("ENST00000233616.8", 2), ("ENST00000448666.6", 6), ("ENST00000452063.7", 6)]
        for t, a in zip(transcripts, answer):
            self.assertEqual((t.id, genome.appris_principality(t)), a)

        # Test use a gene object
        transcripts = genome.appris_transcripts(genome.genes['ENSG00000115275'])
        answer = [("ENST00000233616.8", 2), ("ENST00000448666.6", 6), ("ENST00000452063.7", 6)]
        for t, a in zip(transcripts, answer):
            self.assertEqual((t.id, genome.appris_principality(t)), a)

        self.assertEqual(len(genome.appris_transcripts('ENSG00000239779')), 3)

        self.assertEqual(genome.appris_principality('ENST00000233331', as_string=True), 'PRINCIPAL:1')
        self.assertEqual(genome.appris_principality('ENST00000233615', as_string=True), 'PRINCIPAL:2')
        self.assertEqual(genome.appris_principality('ENST00000233616', as_string=True), 'PRINCIPAL:3')
        self.assertEqual(genome.appris_principality('ENST00000452361', as_string=True), 'PRINCIPAL:5')
        self.assertEqual(genome.appris_principality('ENST00000409737', as_string=True), 'ALTERNATIVE:1')
        self.assertEqual(genome.appris_principality('ENST00000452063', as_string=True), 'ALTERNATIVE:2')

    def test_appris_principality_ucsc_refseq(self):
        genome = MiniGenome('ucsc_refseq.2017-06-25')
        answers = (
            ("NR_037849", None),
            ("NM_012477", 0),
            ("NM_001146158", 6),
            ("NM_006302", 2),
            ("NM_031288", 0),
            ("NM_177552", 0),
            ("NM_001017390", 0),
        )
        for t, p in answers:
            self.assertEqual(genome.appris_principality(t), p)

        gid = '23559'  # NM_012477
        for t in genome.appris_transcripts(gid):
            self.assertTrue(t in genome.genes[gid].transcripts)

        # for ambiguously mapped genes/transcripts, the appris mappings
        # should still be constrained by their genomic range (ie, transcripts
        # are not duplicated across all gene remaps)
        sult1a3 = [x for x in genome.genes if x.name == "SULT1A3"]
        assert len(sult1a3) > 1
        for gene in sult1a3:
            for t in genome.appris_transcripts(gene):
                self.assertTrue(t in gene.transcripts)

        # check sorted by principality
        # SEZ6L2  26470   NM_001243332.1  -   ALTERNATIVE:2
        # SEZ6L2  26470   NM_012410.3 CCDS10658.1 ALTERNATIVE:2
        # SEZ6L2  26470   NM_001243333.1  CCDS58447.1 ALTERNATIVE:2
        # SEZ6L2  26470   NM_001114099.2  -   ALTERNATIVE:2
        # SEZ6L2  26470   NM_201575.3 CCDS10659.1 PRINCIPAL:4
        principalities = [genome.appris_principality(x) for x in genome.appris_transcripts('26470')]
        self.assertEqual(sorted(principalities), principalities)

    def test_mane_simple(self):
        genome = MiniGenome("gencode.v41")
        self.assertTrue(len(genome.mane_transcripts()) > 0)
        self.assertTrue(
            set(genome.mane_transcripts()).issubset(set(genome.transcripts))
        )

    def test_mane_transcripts_unavailable(self):
        genome = MiniGenome("gencode.v29lift37")
        with self.assertRaisesRegex(ValueError, "MANE not supported for annotation"):
            genome.mane_transcripts()

    def test_mane_gencode(self):
        answers = [
            ("ENSG00000115274", "ENST00000233331.12"),
            ("ENSG00000239779", "ENST00000233615.7"),
            ("ENSG00000115275", "ENST00000448666.7"),
        ]

        genome = MiniGenome("gencode.v41")
        for g, t in answers:
            self.assertEqual(genome.mane_transcripts(genome.genes[g]), [genome.transcripts[t]])

        self.assertEqual(
            set(genome.mane_transcripts()),
            {
                genome.transcripts["ENST00000233331.12"],
                genome.transcripts["ENST00000233615.7"],
                genome.transcripts["ENST00000448666.7"],
            },
        )


    def test_mane_ncbi(self):
        answers = [
            ("INO80B", "NM_031288.4"),
            ("WBP1", "NM_012477.4"),
            ("MOGS", "NM_006302.3"),
        ]

        genome = MiniGenome("ncbi_refseq.hg38.p14_RS_2024_08")
        from pprint import pprint
        for g, t in answers:
            gene = genome.genes.first_by_name(g)
            self.assertEqual(genome.mane_transcripts(gene), [genome.transcripts[t]])

        self.assertEqual(
            set(genome.mane_transcripts()),
            {
                genome.transcripts["NM_031288.4"],
                genome.transcripts["NM_012477.4"],
                genome.transcripts["NM_006302.3"],
            },
        )

    def test_cache(self):
        self.assertIs(Genome('hg19'), Genome('hg19'))
        self.assertIsNot(Genome('hg19'), MiniGenome('hg19'))
        self.assertIsNot(MiniGenome('hg19'), MiniGenome('hg19'))

    def test_serialize(self):
        genomes = [
            MiniGenome("gencode.v29"),
            MiniGenome("gencode.v29lift37"),
            MiniGenome("hg19"),
            MiniGenome("hg19.p13.plusMT"),
            MiniGenome("hg38.p12"),
        ]

        for genome in genomes:
            eq = [x for x in genomes if x == genome]
            neq = [x for x in genomes if x != genome]
            self.assertEqual(len(eq), 1)
            self.assertEqual(len(neq), len(genomes) - 1)

        # Check that pickling / unpickling works
        import pickle
        stored = pickle.dumps(genomes, pickle.HIGHEST_PROTOCOL)
        loaded = pickle.loads(stored)

        self.assertListEqual(loaded, genomes)

        g1 = genomes[0]
        g2 = loaded[0]

        self.assertEqual(g1._chromosome_sizes, g2._chromosome_sizes)
        self.assertEqual(g1._appris_indices, g2._appris_indices)
        self.assertEqual(g1._appris_transcripts, g2._appris_transcripts)
        self.assertEqual(g1._appris_transcripts_by_gene, g2._appris_transcripts_by_gene)
        self.assertEqual(g1._appris_principality_strings, g2._appris_principality_strings)


########################################################

if __name__ == "__main__":
    unittest.main()
