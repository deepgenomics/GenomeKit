# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import unittest
import gc
from genome_kit import GenomeDNA
from genome_kit import Interval
from . import MiniGenome

########################################################


class AttrTest(unittest.TestCase):
    def test(self):
        dna = MiniGenome().dna  # Make sure
        self.assertIsInstance(dna, GenomeDNA)
        self.assertIsInstance(dna.filename, str)
        self.assertTrue(dna.filename.endswith(".2bit"))
        with self.assertRaises(TypeError):
            dna.filename = "hg15.2bit"
        with self.assertRaises(AttributeError):
            dna.my_attr = "foo"
        del dna
        gc.collect()
        self.assertEqual(gc.garbage, [])

    def test_arg_typecheck(self):
        genome = MiniGenome()
        with self.assertRaises(TypeError):
            genome.dna(100)


########################################################


class ExtractTest(unittest.TestCase):
    def test(self):
        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        genome = MiniGenome()
        dna = genome.dna

        # Extract the last exon of MOGS transcript ENST00000414701.1
        i = Interval("chr2", "+", 8023, 8039, genome)
        self.assertEqual(dna(i), "CCAGAAGACATTGTAG")
        self.assertEqual(dna(i.as_opposite_strand()), "CTACAATGTCTTCTGG")
        del i

        self.assertEqual(dna(Interval("chr2", '+', 8023, 8023, genome)), "")
        self.assertEqual(dna(Interval("chr2", '+', 8023, 8024, genome)), "C")
        self.assertEqual(dna(Interval("chr2", '+', 8024, 8025, genome)), "C")
        self.assertEqual(dna(Interval("chr2", '+', 8025, 8026, genome)), "A")
        self.assertEqual(dna(Interval("chr2", '+', 8026, 8027, genome)), "G")

        # TODO: test extraction of sequences of all possible lengths
        #       0..100 and at all possible offsets modulo 0..63
        #       to ensure decoding works properly in all cases.

        # Make sure we can't call it with junk arguments
        with self.assertRaises(TypeError):
            dna()

        # Make sure we can't extract coord/interval for wrong reference genome
        with self.assertRaises(ValueError):
            dna(Interval("chr2", "+", 8023, 8023, "hg38.p12"))  # empty interval
        with self.assertRaises(ValueError):
            dna(Interval("chr2", "+", 8023, 8039, "hg38.p12"))

        # Query non-existent chromosome, empty and non-empty intervals
        with self.assertRaises(ValueError):
            dna(Interval("chrX", "+", 0, 0, genome))
        with self.assertRaises(ValueError):
            dna(Interval("chrX", "+", 0, 10, genome))

        del dna

        # Also try a tiny 2bit file containing just two chromosomes
        genome = MiniGenome('test_genome_alt')
        dna = genome.dna
        self.assertEqual(dna(Interval("chr1", "+", 0, 4, genome)), "AAAA")
        self.assertEqual(dna(Interval("chr2", "+", 0, 4, genome)), "CCCC")
        self.assertEqual(dna(Interval("chr2", "-", 0, 4, genome)), "GGGG")
        with self.assertRaises(IndexError):
            dna(Interval("chr1", "+", 0, 5, genome))
        with self.assertRaises(IndexError):
            dna(Interval("chr1", "+", -1, 4, genome))

        gc.collect()
        self.assertEqual(gc.garbage, [])

    def test_allow_outside_chromosome(self):
        genome = MiniGenome()
        dna = genome.dna

        chr2size = genome.chromosome_size("chr2")

        # ranges completely outside of chromosome
        with self.assertRaises(IndexError):
            dna(Interval("chr2", "+", -10, 0, genome), allow_outside_chromosome=True)
        with self.assertRaises(IndexError):
            dna(Interval("chr2", "+", -10, -1, genome), allow_outside_chromosome=True)
        with self.assertRaises(IndexError):
            dna(Interval("chr2", "+", chr2size, chr2size+10, genome), allow_outside_chromosome=True)

        # allow_outside_chromosome=False
        with self.assertRaises(IndexError):
            dna(Interval("chr2", "+", -1, 10, genome), allow_outside_chromosome=False)
        with self.assertRaises(IndexError):
            dna(Interval("chr2", "+", chr2size-10, chr2size+1, genome), allow_outside_chromosome=False)

        # intervals that cross chrom start boundaries
        first50bases = "CGCGAAGCGTGCCCCGCACAAGGATGGTTGCCATGAACCGGAAGTAACTG"
        self.assertEqual(dna(Interval("chr2", "+", 0, 50, genome)), first50bases)
        # test word boundaries (word = 32 bits = 16 bases)
        for j in range(1, 10):
            for i in range(1, 50):
                self.assertEqual(dna(Interval("chr2", "+", -j, i, genome), allow_outside_chromosome=True),
                                 f"{j*'N'}{first50bases[:i]}", f"failed for i={i} j={j}")

        # intervals that cross chrom end boundaries
        last50bases = "GTTAGGTTAAGTCCCAGAGTCCGGAACCGCTGCCTGCGGCTTGACACAGG"
        self.assertEqual(dna(Interval("chr2", "+", chr2size-50, chr2size, genome)), last50bases)
        # test word boundaries
        for j in range(1, 10):
            for i in range(2, 50):
                self.assertEqual(dna(Interval("chr2", "+", chr2size-i, chr2size+j, genome), allow_outside_chromosome=True),
                                 f"{last50bases[-i:]}{j*'N'}", f"failed for i={i} j={j}")

        # cross both start and end boundaries
        actual = dna(Interval("chr2", "+", -10, chr2size+10, genome), allow_outside_chromosome=True)
        self.assertEqual(len(actual), 10 + chr2size + 10)
        self.assertEqual(f"{'N'*10}{first50bases}", actual[:60])
        self.assertEqual(f"{last50bases}{'N'*10}", actual[-60:])



########################################################

if __name__ == "__main__":
    unittest.main()
