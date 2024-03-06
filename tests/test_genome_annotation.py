# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import unittest
import gc
import os
from genome_kit import Genome
from genome_kit import GenomeAnnotation
from genome_kit import GeneTable
from genome_kit import TranscriptTable
from genome_kit import ExonTable
from genome_kit import IntronTable
from genome_kit import CdsTable
from genome_kit import UtrTable
from genome_kit import Interval
from genome_kit import Gene
from genome_kit import Transcript
from genome_kit import Exon
from genome_kit import Intron
from genome_kit import Cds
from genome_kit import Utr
from . import MiniGenome
from . import check_pythonic_indexing


def _check_interval_query(self, elems, query, expected):
    """Given query = (chrom, strand, start, end) and
    expected = (num_end5_within, num_end3_within, num_within, num_overlap)
    runs each type of interval query and checks that the
    results match the expected result from filtering 'elems'
    in a brute force manner.

    Here 'self' is the unittest.TestCase object, and
    num_end5_within, num_end3_within, num_within, num_overlap
    are the expected number of results for each type of query.
    """

    # Brute force filtering functions applied to each element in 'elems'.
    # Each function f(interval) here returns a new function g(elem) that, when called
    # tests whether 'elem' is is within/overlaps the original interval.
    def end5_within(interval):
        return lambda elem: elem.end5.within(interval)

    def end3_within(interval):
        return lambda elem: elem.end3.within(interval)

    def within(interval):
        return lambda elem: elem.within(interval)

    def overlaps(interval):
        return lambda elem: elem.overlaps(interval)

    def exact(interval):
        return lambda elem: elem.interval == interval

    # First assert that the elements can be put in a set and compared
    self.assertEqual(set(elems), set(elems))

    # Convert query to interval
    chrom, strand, start, end = query
    interval = Interval(chrom, strand, start, end, elems[0].refg)

    result = (
        elems.find_5p_within(interval),
        elems.find_3p_within(interval),
        elems.find_within(interval),
        elems.find_overlapping(interval),
        elems.find_exact(interval),
    )
    bruteforce = (
        list(filter(end5_within(interval), elems)),
        list(filter(end3_within(interval), elems)),
        list(filter(within(interval), elems)),
        list(filter(overlaps(interval), elems)),
        list(filter(exact(interval), elems)),
    )
    self.assertTupleEqual(tuple(len(x) for x in result), expected)
    self.assertTupleEqual(tuple(len(x) for x in result), tuple(len(x) for x in bruteforce))
    for i in range(len(result)):
        self.assertSetEqual(set(result[i]), set(bruteforce[i]))


def _check_walk(self, elems, next_attr, prev_attr):
    n = len(elems)
    if n == 0:
        return

    # Test forward iteration
    curr = elems[0]
    for i in range(n):
        self.assertEqual(curr, elems[i])
        curr = getattr(curr, next_attr)
    self.assertIs(curr, None)

    # Test reverse iteration
    curr = elems[-1]
    for i in range(n):
        self.assertEqual(curr, elems[n - i - 1])
        curr = getattr(curr, prev_attr)
    self.assertIs(curr, None)


def _check_updn_stream(self, elems):
    for i in range(len(elems) - 1):
        self.assertTrue(elems[i].upstream_of(elems[i + 1]))
        self.assertTrue(elems[i + 1].dnstream_of(elems[i]))


########################################################


class TestNewDelete(unittest.TestCase):
    def test(self):
        g = MiniGenome()
        a = g.annotation
        del a
        del g
        gc.collect()
        self.assertEqual(gc.garbage, [])

        a = MiniGenome().annotation  # Try holding on to GenomeAnno ref only, make sure genome refcounting works
        del a
        gc.collect()
        self.assertEqual(gc.garbage, [])


class TestCommon(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        cls.genome = MiniGenome("gencode.v29lift37")
        cls.other_genome = MiniGenome("ucsc_refseq.2017-06-25")

    @classmethod
    def tearDownClass(cls):
        del cls.genome

    def test_basic_attributes(self):
        # Check GenomeAnnotation attributes
        self.assertIsInstance(self.genome.annotation, GenomeAnnotation)
        self.assertIsInstance(self.genome.anno, GenomeAnnotation)  # shortname
        self.assertIsInstance(self.genome.anno.genes, GeneTable)
        self.assertIsInstance(self.genome.anno.transcripts, TranscriptTable)
        self.assertIsInstance(self.genome.anno.trans, TranscriptTable)  # shortname
        self.assertIsInstance(self.genome.anno.exons, ExonTable)
        self.assertIsInstance(self.genome.anno.introns, IntronTable)
        self.assertIsInstance(self.genome.anno.cdss, CdsTable)
        self.assertIsInstance(self.genome.anno.utr5s, UtrTable)
        self.assertIsInstance(self.genome.anno.utr3s, UtrTable)
        self.assertIsInstance(self.genome.anno.filename, str)
        self.assertIsInstance(self.genome.anno.__repr__(), str)
        self.assertIsInstance(GenomeAnnotation.binary_version(), int)

        # Check that convenient shorthand attributes on Genome object work
        self.assertIs(self.genome.genes, self.genome.annotation.genes)
        self.assertIs(self.genome.transcripts, self.genome.annotation.transcripts)
        self.assertIs(self.genome.trans, self.genome.annotation.transcripts)
        self.assertIs(self.genome.exons, self.genome.annotation.exons)
        self.assertIs(self.genome.introns, self.genome.annotation.introns)
        self.assertIs(self.genome.cdss, self.genome.annotation.cdss)
        self.assertIs(self.genome.utr5s, self.genome.annotation.utr5s)
        self.assertIs(self.genome.utr3s, self.genome.annotation.utr3s)

        # Check that attributes can't be added dynamically
        with self.assertRaises(TypeError):
            self.genome.annotation.filename = "mygenome.gff3"
        with self.assertRaises(AttributeError):
            self.genome.annotation.my_attr = "foo"
        self.assertFalse(hasattr(self.genome.annotation, "my_attr"))

        # Check invalid gene ID
        with self.assertRaises(KeyError):
            self.genome.genes["ENSG00000239779.6999"]

        # Manually select for the same interval
        gene = self.genome.genes[0]
        self.assertEqual(gene.interval, Interval(gene.chrom, gene.strand, gene.start, gene.end, gene.refg))

        self.assertEqual(self.genome.genes.first_by_name("INO80B"), gene)
        self.assertIsNone(self.genome.genes.first_by_name("doesn't exist"))

        # Check invalid attributes
        with self.assertRaises(AttributeError):
            self.genome.genes.invalid_attr
        with self.assertRaises(AttributeError):
            self.genome.genes[0].invalid_attr

        # Check adding/removing a dynamic attribute
        with self.assertRaises(AttributeError):
            self.genome.genes.invalid_attr = "hello"
        with self.assertRaises(AttributeError):
            self.genome.genes[0].invalid_attr = "hello"
        self.assertFalse(hasattr(self.genome.genes, "invalid_attr"))
        self.assertFalse(hasattr(self.genome.genes[0], "invalid_attr"))

        # Check __init__ doesn't crash on internally derived intervals
        self.check_invalid_init(Gene)
        self.check_invalid_init(Transcript)
        self.check_invalid_init(Exon)
        self.check_invalid_init(Intron)
        self.check_invalid_init(Cds)
        self.check_invalid_init(Utr)

    def check_invalid_init(self, type):
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            type('chr1', '+', 0, 0, self.genome)

    def test_pythonic_indexing(self):

        # Check that genes table can be indexed by integer in the usual Python-esque way
        check_pythonic_indexing(self, self.genome.genes)
        check_pythonic_indexing(self, self.genome.transcripts)
        check_pythonic_indexing(self, self.genome.exons)
        check_pythonic_indexing(self, self.genome.introns)
        check_pythonic_indexing(self, self.genome.cdss)
        check_pythonic_indexing(self, self.genome.utr5s)
        check_pythonic_indexing(self, self.genome.utr3s)

    def test_element_equality(self):
        def check(items):
            # Each item should be equal only to itself, and != with all others
            for item in items:
                eq = [x for x in items if x == item]
                neq = [x for x in items if x != item]
                self.assertEqual(len(eq), 1)
                self.assertEqual(len(neq), len(items) - 1)

        check(self.genome.genes)
        check(self.genome.transcripts)
        check(self.genome.exons)
        check(self.genome.introns)
        check(self.genome.cdss)
        check(self.genome.utr5s)
        check(self.genome.utr3s)

        self.assertNotEqual(self.genome.trans[0].exons[0], self.other_genome.trans[0].exons[0])

    def test_element_hash(self):
        def check(items):
            # Each item should have a unique hash, with the current hash functions
            # and the current MiniGenome. If the genome or hash functions change
            # in future, then there may be legitimate hash collisions.
            for item in items:
                eq = [x for x in items if hash(x) == hash(item)]
                neq = [x for x in items if hash(x) != hash(item)]
                self.assertEqual(len(eq), 1)
                self.assertEqual(len(neq), len(items) - 1)

        check(self.genome.genes)
        check(self.genome.transcripts)
        check(self.genome.exons)
        check(self.genome.introns)
        check(self.genome.cdss)
        check(self.genome.utr5s)
        check(self.genome.utr3s)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full dgannos on CI")
    def test_serialize(self):
        genome = Genome("ucsc_refseq.2017-06-25")
        tables = [
            getattr(genome, table)
            for table in ["genes", "trans", "exons", "introns", "cdss", "utr5s", "utr3s"]
        ]
        features_by_type = [
            [
                x[0],
                x[-1],
                x[len(x) // 2],
            ]
            for x in tables
        ]

        # Check that pickling / unpickling works
        import pickle
        stored = pickle.dumps(features_by_type, pickle.HIGHEST_PROTOCOL)
        loaded = pickle.loads(stored)
        self.assertEqual(loaded, features_by_type)

########################################################


class TestGencode(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        cls.genome = MiniGenome("gencode.v29lift37")

    @classmethod
    def tearDownClass(cls):
        del cls.genome

    def test_gene_interval_queries(self):
        def check(query, expected):
            _check_interval_query(self, self.genome.genes, query, expected)

        # yapf: disable
        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        # expected = (num_end5_within, num_end3_within, num_within, num_overlap, num_exact)
        check(("chr2", "-",     0, 11000), (1, 1, 1, 1, 0))  # - strand (MOGS)
        check(("chr2", "+",     0, 11000), (3, 3, 3, 3, 0))  # + strand (INO80B, INO80B-WBP1, WBP1)
        check(("chr2", "+",    49,  2987), (2, 1, 1, 2, 1))  # INO80B exactly
        check(("chr2", "+",    50,  2987), (1, 1, 0, 2, 0))  # INO80B +1 5p end
        check(("chr2", "+",    49,  2986), (2, 0, 0, 2, 0))  # INO80B -1 3p end
        check(("chr2", "-",  6083, 10437), (1, 1, 1, 1, 1))  # MOGS exactly
        check(("chr2", "-",  6083, 10436), (0, 1, 0, 1, 0))  # MOGS -1 5p end
        check(("chr2", "-",  6084, 10437), (1, 0, 0, 1, 0))  # MOGS +1 3p end
        # yapf: enable

    def test_transcript_interval_queries(self):
        def check(query, expected):
            _check_interval_query(self, self.genome.transcripts, query, expected)

        # yapf: disable
        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        # expected = (num_end5_within, num_end3_within, num_within, num_overlap, num_exact)
        check(("chr2", "-",     0, 11000), (20, 20, 20, 20, 0))  # - strand (MOGS)
        check(("chr2", "+",     0, 11000), (25, 25, 25, 25, 0))  # + strand (INO80B, INO80B-WBP1, WBP1)
        check(("chr2", "+",  3495,  5883), ( 7, 12,  7, 16, 1))  # ENST00000409737.5 exactly
        check(("chr2", "+",  3496,  5883), ( 6, 12,  6, 16, 0))  # ENST00000409737.5 +1 5p end
        check(("chr2", "+",  3495,  5882), ( 7, 11,  6, 16, 0))  # ENST00000409737.5 -1 3p end
        check(("chr2", "-",  6083, 10409), (10, 20, 10, 20, 1))  # ENST00000448666.6 exactly
        check(("chr2", "-",  6083, 10408), ( 9, 20,  9, 20, 0))  # ENST00000448666.6 -1 5p end
        check(("chr2", "-",  6084, 10409), (10, 18,  9, 20, 0))  # ENST00000448666.6 +1 3p end
        # yapf: enable

    def test_exon_interval_queries(self):
        def check(query, expected):
            _check_interval_query(self, self.genome.exons, query, expected)

        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        # expected = (num_5p_within, num_3p_within, num_within, num_overlap, num_exact)
        n = len([e for e in self.genome.exons if e.chrom == 'chr2' and e.strand == '-'])
        p = len([e for e in self.genome.exons if e.chrom == 'chr2' and e.strand == '+'])
        check(("chr2", "-", 0, 11000), (n, n, n, n, 0))  # - strand (MOGS)
        check(("chr2", "+", 0, 11000), (p, p, p, p, 0))  # + strand (INO80B, INO80B-WBP1, WBP1)
        # TODO: test regions around specific exons.

    def test_find_invalid_refg(self):
        trans = self.genome.annotation.transcripts
        interval = Interval('chr1', '+', 0, 0, 'hg38.p12')

        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_3p_aligned(interval)
        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_3p_within(interval)
        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_5p_aligned(interval)
        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_5p_within(interval)
        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_exact(interval)
        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_overlapping(interval)
        with self.assertRaisesRegex(ValueError, "Reference"):
            trans.find_within(interval)

    def test_find_aligned(self):
        anno = self.genome.annotation

        # yapf: disable
        # For mini genome chr2_74682100_74692599_h37, pull out an exon
        # on each strand and make sure find_5/3p_aligned works correctly
        # for both the ExonTable and IntronTable.
        exon = anno.transcripts['ENST00000233331.11'].exons[3]  # + strand exon
        self.assertEqual(len(anno.exons.find_5p_aligned(exon     )), 10)
        self.assertEqual(len(anno.exons.find_5p_aligned(exon.end5)), 10)
        self.assertEqual(len(anno.exons.find_5p_aligned(exon.end3)), 0)
        self.assertEqual(len(anno.exons.find_3p_aligned(exon     )), 11)
        self.assertEqual(len(anno.exons.find_3p_aligned(exon.end5)), 0)
        self.assertEqual(len(anno.exons.find_3p_aligned(exon.end3)), 11)
        self.assertEqual(len(anno.introns.find_5p_aligned(exon     )), 0)
        self.assertEqual(len(anno.introns.find_5p_aligned(exon.end5)), 0)
        self.assertEqual(len(anno.introns.find_5p_aligned(exon.end3)), 11)
        self.assertEqual(len(anno.introns.find_3p_aligned(exon     )), 0)
        self.assertEqual(len(anno.introns.find_3p_aligned(exon.end5)), 10)
        self.assertEqual(len(anno.introns.find_3p_aligned(exon.end3)), 0)

        exon = anno.transcripts['ENST00000233616.8'].exons[2]  # - strand exon
        self.assertEqual(len(anno.exons.find_5p_aligned(exon     )), 13)
        self.assertEqual(len(anno.exons.find_5p_aligned(exon.end5)), 13)
        self.assertEqual(len(anno.exons.find_5p_aligned(exon.end3)), 0)
        self.assertEqual(len(anno.exons.find_3p_aligned(exon     )), 16)
        self.assertEqual(len(anno.exons.find_3p_aligned(exon.end5)), 0)
        self.assertEqual(len(anno.exons.find_3p_aligned(exon.end3)), 16)
        self.assertEqual(len(anno.introns.find_5p_aligned(exon     )), 0)
        self.assertEqual(len(anno.introns.find_5p_aligned(exon.end5)), 0)
        self.assertEqual(len(anno.introns.find_5p_aligned(exon.end3)), 16)
        self.assertEqual(len(anno.introns.find_3p_aligned(exon     )), 0)
        self.assertEqual(len(anno.introns.find_3p_aligned(exon.end5)), 13)
        self.assertEqual(len(anno.introns.find_3p_aligned(exon.end3)), 0)
        # yapf: enable

    def test_gene_attributes(self):
        # Pick a known gene, check all its attributes
        gene = self.genome.genes["ENSG00000115275.12"]  # MOGS gene
        self.assertEqual(gene, self.genome.genes["ENSG00000115275"])  # ID without version
        self.assertEqual(gene.id, "ENSG00000115275.12")
        self.assertEqual(gene.name, "MOGS")
        self.assertEqual(gene.type, "protein_coding")
        self.assertEqual(gene.level, 1)
        self.assertEqual(gene.interval, Interval("chr2", "-", 6083, 10437, self.genome))
        self.assertEqual(gene.trans, gene.transcripts)  # shortname
        self.assertEqual(len(gene.transcripts), 20)
        self.assertIsInstance(gene, Interval)
        self.assertIsInstance(gene.__repr__(), str)

    def test_transcript_attributes(self):
        # Pull out a gene's known transcripts, and check their IDS
        gene = self.genome.genes["ENSG00000115275.12"]  # MOGS gene
        trans = gene.transcripts
        self.assertEqual(trans, [self.genome.transcripts[t.id] for t in trans])
        self.assertEqual(len(trans), 20)
        self.assertListEqual(
            [x.id for x in (trans[1], trans[12], trans[-1])],
            ["ENST00000233616.8", "ENST00000452063.7", "ENST00000486036.1"],
        )

        # Check one of the transcripts' attributes
        tran = trans[1]
        self.assertEqual(tran, self.genome.transcripts['ENST00000233616'])  # Without version
        self.assertEqual(tran.type, "protein_coding")
        self.assertEqual(tran.level, 2)
        self.assertEqual(tran.tsl, 1)
        self.assertEqual(tran.ccds_id, "CCDS42700.1")
        self.assertEqual(tran.protein_id, "ENSP00000233616.4")
        self.assertEqual(tran.interval, Interval("chr2", "-", 6083, 10437, self.genome))
        self.assertEqual(tran.gene, gene)  # Check parent is the original gene
        self.assertIs(tran.product, None)  # Not specified by GENCODE
        self.assertEqual(len(tran.exons), 4)
        self.assertEqual(len(tran.introns), 3)
        self.assertEqual(len(tran.cdss), 4)
        self.assertEqual(len(tran.utr5s), 1)
        self.assertEqual(len(tran.utr3s), 1)
        self.assertIsInstance(tran, Interval)
        self.assertIsInstance(tran.__repr__(), str)

    def test_exon_attributes(self):
        # Pull out the transcript's exons
        gene = self.genome.genes["ENSG00000115275.12"]  # MOGS gene
        tran = gene.trans[1]
        exons = tran.exons
        self.assertListEqual(
            [x.id for x in exons],
            [
                "ENSE00000846737.3",
                "ENSE00003462747.1",
                "ENSE00003542785.1",
                "ENSE00001853789.1",
            ],
        )
        for index, exon in enumerate(exons):
            self.assertEqual(exon.index, index)
            self.assertEqual(exon.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(exon.tran, tran)  # shortname
            self.assertIsInstance(exon, Interval)
            self.assertIsInstance(exon.__repr__(), str)

        # Walk the transcript's exons as a doubly linked list
        _check_walk(self, exons, "next_exon", "prev_exon")
        _check_updn_stream(self, exons)

        # Check a specific exon's attributes in more detail
        exon = exons[0]
        self.assertEqual(exon.interval, Interval("chr2", "-", 9922, 10437, self.genome))

        # Check that CDS of a non-coding exon is None
        tran = self.genome.transcripts["ENST00000452063.7"]
        exon = tran.exons[0]
        self.assertIsNone(exon.cds)
        self.assertEqual(exon.interval, exon.utr5.interval)
        self.assertIsNone(exon.utr3)

    def test_intron_attributes(self):
        # Pull out the transcript's introns
        gene = self.genome.genes["ENSG00000115275.12"]  # MOGS gene
        tran = gene.trans[0]
        introns = tran.introns
        self.assertEqual(len(introns), 3)
        for index, intron in enumerate(introns):
            self.assertEqual(intron.index, index)
            self.assertEqual(intron.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(intron.tran, tran)  # shortname
            self.assertIsInstance(intron, Interval)
            self.assertIsInstance(intron.__repr__(), str)

            # Check that exon.next_intron and intron.next_exon work
            self.assertEqual(intron.prev_exon, tran.exons[index])
            self.assertEqual(intron.next_exon, tran.exons[index + 1])
            self.assertEqual(intron.prev_exon.next_intron, intron)
            self.assertEqual(intron.next_exon.prev_intron, intron)

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, introns, "next_intron", "prev_intron")
        _check_updn_stream(self, introns)

        # Check a specific intron's attributes
        intron = introns[0]
        self.assertEqual(intron.interval, Interval("chr2", "-", 9749, 9922, self.genome))

    def test_cds_attributes(self):
        # Pull out the transcript's CDSs
        gene = self.genome.genes["ENSG00000115275.12"]  # MOGS gene
        tran = gene.trans[0]
        exon = tran.exons[0]
        cdss = tran.cdss
        phase = 0
        self.assertEqual(len(cdss), 4)
        for index, cds in enumerate(cdss):
            self.assertEqual(cds, tran.exons[index].cds)  # This transcript has one CDS per exon
            self.assertEqual(cds.exon.index, index)  # This transcript has one CDS per exon
            self.assertEqual(cds.phase, phase)
            self.assertEqual(cds.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(cds.tran, tran)  # shortname
            self.assertIsInstance(cds, Interval)
            self.assertIsInstance(cds.__repr__(), str)
            phase = 2 - ((2 - phase) + len(cds)) % 3  # phase = num bp to skip before start of next codon

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, cdss, "next_cds", "prev_cds")
        _check_updn_stream(self, cdss)

        # Check a particular CDS's attributes in more detail
        cds = tran.cdss[0]
        self.assertEqual(exon, exon.cds.exon)  # Make sure exon and cds refer to one another properly
        self.assertEqual(cds.phase, 0)
        self.assertEqual(cds.interval, Interval("chr2", "-", 9922, 10274, self.genome))
        self.assertEqual(cds.exon, exon)

    def test_utr5_attributes(self):
        gene = self.genome.genes["ENSG00000149927.17"]  # DOC2A gene
        tran = gene.trans[7]
        utr5s = tran.utr5s
        self.assertEqual(len(utr5s), 4)
        for index, utr5 in enumerate(utr5s):
            self.assertEqual(utr5, tran.exons[index].utr5)
            self.assertEqual(utr5.exon.index, index)
            self.assertEqual(utr5.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(utr5.tran, tran)  # shortname
            self.assertIsInstance(utr5, Interval)
            self.assertIsInstance(utr5.__repr__(), str)

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, utr5s, "next_utr", "prev_utr")
        _check_updn_stream(self, utr5s)

        exon = tran.exons[0]
        utr5 = tran.utr5s[0]
        # Make sure exon and utr5 refer to one another properly
        self.assertEqual(exon, exon.utr5.exon)
        self.assertEqual(utr5.exon, exon)

    def test_utr3_attributes(self):
        gene = self.genome.genes["ENSG00000149927.17"]  # DOC2A gene
        tran = gene.trans[7]
        utr3s = tran.utr3s
        self.assertEqual(len(utr3s), 5)
        for index, utr3 in enumerate(utr3s):
            idx = index + 5 # first 5 exons don't have UTR3s
            self.assertEqual(utr3, tran.exons[idx].utr3)
            self.assertEqual(utr3.exon.index, idx)
            self.assertEqual(utr3.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(utr3.tran, tran)  # shortname
            self.assertIsInstance(utr3, Interval)
            self.assertIsInstance(utr3.__repr__(), str)

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, utr3s, "next_utr", "prev_utr")
        _check_updn_stream(self, utr3s)

        exon = tran.exons[9]
        utr3 = tran.utr3s[4]
        self.assertEqual(exon, exon.utr3.exon)  # Make sure exon and utr3 refer to one another properly
        self.assertEqual(utr3.exon, exon)

    def test_utr5_boundaries(self):
        tran = self.genome.transcripts["ENST00000561671.5"]  # DOC2A gene
        utr5 = tran.utr5s[0]
        # Negative strand
        # chr16	HAVANA	five_prime_UTR	551018	551142	.	-
        self.assertEqual(Interval("chr16", "-", 551017, 551142, self.genome), utr5.interval)
        # Positive strand
        # chr16	HAVANA	five_prime_UTR	352203	352438	.	+
        utr5 = self.genome.transcripts["ENST00000358758.11"].utr5s[0]
        self.assertEqual(Interval("chr16", "+", 352202, 352438, self.genome), utr5.interval)

        # partially coding exon (positive strand)
        tran = self.genome.transcripts["ENST00000358758.11"]
        exon = tran.exons[1]
        self.assertEqual(exon.start, exon.utr5.start)
        self.assertNotEqual(exon.end, exon.utr5.end)
        self.assertEqual(exon.end, exon.cds.end)
        self.assertEqual(exon.utr5.end, exon.cds.start)

        # partially coding exon (negative strand)
        tran = self.genome.transcripts["ENST00000452063.7"]
        exon = tran.exons[1]
        self.assertEqual(exon.end, exon.utr5.end)
        self.assertNotEqual(exon.start, exon.utr5.start)
        self.assertEqual(exon.start, exon.cds.start)
        self.assertEqual(exon.cds.end, exon.utr5.start)

        # non-coding exon
        tran = self.genome.transcripts["ENST00000452063.7"]
        exon = tran.exons[0]
        self.assertIs(exon.cds, None)
        self.assertEqual(exon.interval, exon.utr5.interval)

    def test_utr3_boundaries(self):
        tran = self.genome.transcripts["ENST00000561671.5"]  # DOC2A gene
        utr3 = tran.utr3s[4]

        # Negative strand
        # chr16	HAVANA	three_prime_UTR	546409	546624	.	-
        self.assertEqual(Interval("chr16", "-", 546408, 546624, self.genome), utr3.interval)
        # Positive strand
        # chr16	ENSEMBL	three_prime_UTR	625048	625471	.	+
        utr5 = self.genome.transcripts["ENST00000627746.2"].utr3s[4]
        self.assertEqual(Interval("chr16", "+", 625047, 625471, self.genome), utr5.interval)

        # partially coding exon (positive strand)
        tran = self.genome.transcripts["ENST00000627746.2"]
        exon = tran.exons[3]
        self.assertEqual(exon.start, exon.cds.start)
        self.assertNotEqual(exon.start, exon.utr3.start)
        self.assertEqual(exon.end, exon.utr3.end)
        self.assertEqual(exon.cds.end, exon.utr3.start)

        # partially coding exon (negative strand)
        tran = self.genome.transcripts["ENST00000561671.5"]
        exon = tran.exons[5]
        self.assertEqual(exon.start, exon.utr3.start)
        self.assertNotEqual(exon.end, exon.utr3.end)
        self.assertEqual(exon.end, exon.cds.end)
        self.assertEqual(exon.utr3.end, exon.cds.start)

        # non-coding exon
        tran = self.genome.transcripts["ENST00000627746.2"]
        exon = tran.exons[7]
        self.assertIs(exon.cds, None)
        self.assertEqual(exon.interval, exon.utr3.interval)

    def test_utr_undifferentiated(self):
        # UTRs marked as `UTR` (rather than `five_prime_UTR` or `three_prime_UTR`)
        tran = self.genome.transcripts["ENST00000551448.1"]  # RP11-231C14.4 gene
        self.assertEqual(2, len(tran.utr5s))
        self.assertEqual(7, len(tran.utr3s))
        self.assertEqual(tran.exons[0].interval, tran.utr5s[0].interval)
        self.assertEqual(tran.exons[1].end, tran.utr5s[1].end)
        # no coding-only exons
        self.assertEqual(tran.exons[2].start, tran.utr3s[0].start)
        self.assertEqual(tran.exons[3].interval, tran.utr3s[1].interval)
        self.assertEqual(tran.exons[8].interval, tran.utr3s[6].interval)



########################################################


class TestUCSCRefSeq(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.genome = MiniGenome("ucsc_refseq.2017-06-25")

    @classmethod
    def tearDownClass(cls):
        del cls.genome

    def test_gene_interval_queries(self):
        def check(query, expected):
            _check_interval_query(self, self.genome.genes, query, expected)

        # yapf: disable
        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        # expected = (num_end5_within, num_end3_within, num_within, num_overlap, num_exact)
        check(("chr2", "-",     0, 11000), (1, 1, 1, 1, 0))  # - strand (MOGS)
        check(("chr2", "+",     0, 11000), (3, 3, 3, 3, 0))  # + strand (INO80B, INO80B-WBP1, WBP1)
        check(("chr2", "+",    49,  2987), (2, 1, 1, 2, 1))  # INO80B exactly
        check(("chr2", "+",    50,  2987), (0, 1, 0, 2, 0))  # INO80B +1 5p end
        check(("chr2", "+",    49,  2986), (2, 0, 0, 2, 0))  # INO80B -1 3p end
        check(("chr2", "-",  6083, 10437), (1, 1, 1, 1, 1))  # MOGS exactly
        check(("chr2", "-",  6083, 10436), (0, 1, 0, 1, 0))  # MOGS -1 5p end
        check(("chr2", "-",  6084, 10437), (1, 0, 0, 1, 0))  # MOGS +1 3p end
        # yapf: enable

    def test_transcript_interval_queries(self):
        def check(query, expected):
            _check_interval_query(self, self.genome.transcripts, query, expected)

        # yapf: disable
        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        # expected = (num_end5_within, num_end3_within, num_within, num_overlap, num_exact)
        check(("chr2", "-",     0, 11000), (2, 2, 2, 2, 0))  # - strand (MOGS)
        check(("chr2", "+",     0, 11000), (3, 3, 3, 3, 0))  # + strand (INO80B, INO80B-WBP1, WBP1)
        check(("chr2", "+",  3426,  5918), (1, 2, 1, 2, 1))  # NM_012477 exactly
        check(("chr2", "+",  3427,  5918), (0, 2, 0, 2, 0))  # NM_012477 +1 5p end
        check(("chr2", "+",  3426,  5917), (1, 0, 0, 2, 0))  # NM_012477 -1 3p end
        check(("chr2", "-",  6083, 10437), (2, 2, 2, 2, 2))  # NM_001146158 exactly
        check(("chr2", "-",  6083, 10436), (0, 2, 0, 2, 0))  # NM_001146158 -1 5p end
        check(("chr2", "-",  6084, 10437), (2, 0, 0, 2, 0))  # NM_001146158 +1 3p end
        # yapf: enable

    def test_exon_interval_queries(self):
        def check(query, expected):
            _check_interval_query(self, self.genome.exons, query, expected)

        # For mini genome chr2_74682100_74692599_h37 (INO80B, INO80B-WBP1, WBP1, MOGS).
        # expected = (num_5p_within, num_3p_within, num_within, num_overlap, num_exact)
        n = len([e for e in self.genome.exons if e.chrom == 'chr2' and e.strand == '-'])
        p = len([e for e in self.genome.exons if e.chrom == 'chr2' and e.strand == '+'])
        check(("chr2", "-", 0, 11000), (n, n, n, n, 0))  # - strand (MOGS)
        check(("chr2", "+", 0, 11000), (p, p, p, p, 0))  # + strand (INO80B, INO80B-WBP1, WBP1)
        # TODO: test regions around specific exons.

    def test_gene_attributes(self):
        # Pick a known gene, check all its attributes
        gene = self.genome.genes["7841"]  # MOGS gene Entrez ID
        self.assertEqual(gene.id, "7841")
        self.assertEqual(gene.name, "MOGS")
        self.assertEqual(gene.type, "protein_coding")
        self.assertEqual(gene.level, None)
        self.assertEqual(gene.interval, Interval("chr2", "-", 6083, 10437, self.genome))
        self.assertEqual(gene.trans, gene.transcripts)  # shortname
        self.assertEqual(len(gene.transcripts), 2)
        self.assertIsInstance(gene, Interval)
        self.assertIsInstance(gene.__repr__(), str)

        # Manually select for the same interval
        self.assertEqual(gene.interval, Interval(gene.chrom, gene.strand, gene.start, gene.end, gene.refg))

    def test_transcript_attributes(self):
        # Pull out a gene's known transcripts, and check their IDS
        gene = self.genome.genes["7841"]  # MOGS gene
        trans = gene.transcripts
        self.assertEqual(trans, [self.genome.transcripts[t.id] for t in trans])
        self.assertEqual(len(trans), 2)
        self.assertEqual(trans[0].id, "NM_001146158")
        self.assertEqual(trans[1].id, "NM_006302")

        # Check one of the transcripts' attributes
        tran = trans[1]
        self.assertEqual(tran.type, "protein_coding")
        self.assertEqual(tran.level, None)
        self.assertEqual(tran.tsl, None)
        self.assertEqual(tran.ccds_id, None)
        self.assertEqual(tran.protein_id, "NP_006293")
        self.assertEqual(tran.interval, Interval("chr2", "-", 6083, 10437, self.genome))
        self.assertEqual(tran.gene, gene)  # Check parent is the original gene
        self.assertEqual(tran.product, "mannosyl-oligosaccharide glucosidase isoform 1")
        self.assertEqual(len(tran.exons), 4)
        self.assertEqual(len(tran.introns), 3)
        self.assertEqual(len(tran.cdss), 4)
        self.assertEqual(len(tran.utr5s), 1)
        self.assertEqual(len(tran.utr3s), 1)
        self.assertIsInstance(tran, Interval)
        self.assertIsInstance(tran.__repr__(), str)

    def test_exon_attributes(self):
        # Pull out the transcript's exons
        gene = self.genome.genes["7841"]  # MOGS gene
        tran = gene.trans[1]
        exons = tran.exons
        self.assertEqual(len(exons), 4)
        for index, exon in enumerate(exons):
            self.assertEqual(exon.id, None)
            self.assertEqual(exon.index, index)
            self.assertEqual(exon.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(exon.tran, tran)  # shortname
            self.assertIsInstance(exon, Interval)
            self.assertIsInstance(exon.__repr__(), str)

        # Walk the transcript's exons as a doubly linked list
        _check_walk(self, exons, "next_exon", "prev_exon")
        _check_updn_stream(self, exons)

        # Check a specific exon's attributes in more detail
        exon = exons[0]
        self.assertEqual(exon.interval, Interval("chr2", "-", 9922, 10437, self.genome))

        # Check that CDS of a non-coding exon is None
        tran = self.genome.transcripts["NM_001146158"]
        exon = tran.exons[0]
        self.assertIsNone(exon.cds)
        self.assertIsNotNone(exon.utr5)
        self.assertIsNone(exon.utr3)

        # non-coding transcript
        self.assertEqual([x for x in self.genome.transcripts['NR_037849'].exons if x.cds], [])
        self.assertEqual([x for x in self.genome.transcripts['NR_037849'].exons if x.utr5], [])
        self.assertEqual([x for x in self.genome.transcripts['NR_037849'].exons if x.utr3], [])

    def test_intron_attributes(self):
        # Pull out the transcript's introns
        gene = self.genome.genes["7841"]  # MOGS gene
        tran = gene.trans[1]
        introns = tran.introns
        self.assertEqual(len(introns), 3)
        for index, intron in enumerate(introns):
            self.assertEqual(intron.index, index)
            self.assertEqual(intron.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(intron.tran, tran)  # shortname
            self.assertIsInstance(intron, Interval)
            self.assertIsInstance(intron.__repr__(), str)

            # Check that exon.next_intron and intron.next_exon work
            self.assertEqual(intron.prev_exon, tran.exons[index])
            self.assertEqual(intron.next_exon, tran.exons[index + 1])
            self.assertEqual(intron.prev_exon.next_intron, intron)
            self.assertEqual(intron.next_exon.prev_intron, intron)

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, introns, "next_intron", "prev_intron")
        _check_updn_stream(self, introns)

        # Check a specific intron's attributes
        intron = introns[0]
        self.assertEqual(intron.interval, Interval("chr2", "-", 9749, 9922, self.genome))

    def test_cds_attributes(self):
        # Pull out the transcript's CDSs
        gene = self.genome.genes["7841"]  # MOGS gene
        tran = gene.trans[1]
        exon = tran.exons[0]
        cdss = tran.cdss
        phase = 0
        self.assertEqual(len(cdss), 4)
        for index, cds in enumerate(cdss):
            self.assertEqual(cds, tran.exons[index].cds)  # This transcript has one CDS per exon
            self.assertEqual(cds.exon.index, index)  # This transcript has one CDS per exon
            self.assertEqual(cds.phase, phase)
            self.assertEqual(cds.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(cds.tran, tran)  # shortname
            self.assertIsInstance(cds, Interval)
            self.assertIsInstance(cds.__repr__(), str)
            phase = 2 - ((2 - phase) + len(cds)) % 3  # phase = num bp to skip before start of next codon

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, cdss, "next_cds", "prev_cds")
        _check_updn_stream(self, cdss)

        # Check a particular CDS's attributes in more detail
        cds = tran.cdss[0]
        self.assertEqual(exon, exon.cds.exon)  # Make sure exon and cds refer to one another properly
        self.assertEqual(cds.phase, 0)
        self.assertEqual(cds.interval, Interval("chr2", "-", 9922, 10274, self.genome))
        self.assertEqual(cds.exon, exon)

        self.assertFalse(self.genome.transcripts['NR_037849'].cdss)

    def test_utr5_attributes(self):
        # Pull out the transcript's UTR5s
        gene = self.genome.genes["124446"]  # TMEM219 gene
        tran = gene.trans[1]
        exon = tran.exons[0]
        utr5s = tran.utr5s
        self.assertEqual(len(utr5s), 2)
        for index, utr5 in enumerate(utr5s):
            self.assertEqual(utr5, tran.exons[index].utr5)  # This transcript has one CDS per exon
            self.assertEqual(utr5.exon.index, index)  # This transcript has one CDS per exon
            self.assertEqual(utr5.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(utr5.tran, tran)  # shortname
            self.assertIsInstance(utr5, Interval)
            self.assertIsInstance(utr5.__repr__(), str)

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, utr5s, "next_utr", "prev_utr")
        _check_updn_stream(self, utr5s)

        # Check a particular UTR5's attributes in more detail
        utr5 = tran.utr5s[0]
        self.assertEqual(exon, exon.utr5.exon)  # Make sure exon and utr5 refer to one another properly
        self.assertEqual(Interval("chr16", "+", 502144, 502247, self.genome), utr5.interval)
        self.assertEqual(utr5.exon, exon)

        # non-coding transcript
        self.assertFalse(self.genome.transcripts['NR_037849'].utr5s)

    def test_utr5_boundaries(self):
        tran = self.genome.transcripts["NM_001286586"]  # CDIPT gene
        utr5 = tran.utr5s[0]
        # Negative strand
        # NM_001286586	chr16	-	398470	403403	399303	401357	6	398470,399556,400695,401220,402700,402929,	399449,399638,400777,401374,402750,403403,
        self.assertEqual(Interval("chr16", "-", 402929, 403403, self.genome), utr5.interval)
        # Positive strand
        # NM_001083613	chr16	+	502144	513167	503258	511660	6	502144,503221,503514,508139,511522,513064,	502247,503423,503704,508369,511693,513167,
        utr5 = self.genome.transcripts["NM_001083613"].utr5s[0] # TMEM219 gene
        self.assertEqual(Interval("chr16", "+", 502144, 502247, self.genome), utr5.interval)

        # partially coding exon (positive strand)
        tran = self.genome.transcripts["NM_001083613"]
        exon = tran.exons[1]
        self.assertEqual(exon.start, exon.utr5.start)
        self.assertNotEqual(exon.end, exon.utr5.end)
        self.assertEqual(exon.end, exon.cds.end)
        self.assertEqual(exon.utr5.end, exon.cds.start)

        # partially coding exon (negative strand)
        tran = self.genome.transcripts["NM_001286586"]
        exon = tran.exons[2]
        self.assertEqual(exon.end, exon.utr5.end)
        self.assertNotEqual(exon.start, exon.utr5.start)
        self.assertEqual(exon.start, exon.cds.start)
        self.assertEqual(exon.cds.end, exon.utr5.start)

        # non-coding exon
        tran = self.genome.transcripts["NM_001083613"]
        exon = tran.exons[0]
        self.assertIs(exon.cds, None)
        self.assertEqual(exon.interval, exon.utr5.interval)

        # UTR ending on exon boundary
        tran = self.genome.transcripts["NM_175900"]
        self.assertEqual(1, len(tran.utr5s))
        self.assertEqual(tran.exons[0].utr5.interval, tran.exons[0].interval)
        self.assertIsNone(tran.exons[0].cds)
        self.assertIsNone(tran.exons[1].utr5)

        # exon that includes UTR5 and UTR3
        tran = self.genome.transcripts["NM_001256443"]
        self.assertIsNotNone(tran.exons[1].utr5)
        self.assertGreater(len(tran.exons[1].utr5), 0)
        self.assertIsNotNone(tran.exons[1].cds)
        self.assertGreater(len(tran.exons[1].cds), 0)
        self.assertIsNotNone(tran.exons[1].utr3)
        self.assertGreater(len(tran.exons[1].utr3), 0)

    def test_utr3_attributes(self):
        # Pull out the transcript's UTR3s
        gene = self.genome.genes["124446"]  # TMEM219 gene
        tran = gene.trans[1]
        utr3s = tran.utr3s
        self.assertEqual(len(utr3s), 2)
        for index, utr3 in enumerate(utr3s):
            idx = index + 4 # first 4 exons don't have UTR3s
            self.assertEqual(utr3, tran.exons[idx].utr3)  # This transcript has one CDS per exon
            self.assertEqual(utr3.exon.index, idx)  # This transcript has one CDS per exon
            self.assertEqual(utr3.transcript, tran)  # Check parent is the original transcript
            self.assertEqual(utr3.tran, tran)  # shortname
            self.assertIsInstance(utr3, Interval)
            self.assertIsInstance(utr3.__repr__(), str)

        # Walk a transcript's introns as a doubly linked list
        _check_walk(self, utr3s, "next_utr", "prev_utr")
        _check_updn_stream(self, utr3s)

        # Check a particular UTR3's attributes in more detail
        utr3 = tran.utr3s[1]
        exon = tran.exons[5]
        self.assertEqual(exon, exon.utr3.exon)  # Make sure exon and utr3 refer to one another properly
        self.assertEqual(Interval("chr16", "+", 513064, 513167, self.genome), utr3.interval)
        self.assertEqual(utr3.exon, exon)

        # non-coding transcript
        self.assertFalse(self.genome.transcripts['NR_037849'].utr3s)

    def test_utr3_boundaries(self):
        tran = self.genome.transcripts["NM_001109891"]  # MAPK3 gene
        utr3 = tran.utr3s[1]
        # Negative strand
        # NM_001109891	chr16	-	654219	663424	656782	663324	8	654219,656750,657008,657784,658161,658463,661938,663154,	654823,656905,657118,657899,658278,658653,662121,663424,
        self.assertEqual(Interval("chr16", "-", 654219, 654823, self.genome), utr3.interval)
        # Positive strand
        # NM_001083613	chr16	+	502144	513167	503258	511660	6	502144,503221,503514,508139,511522,513064,	502247,503423,503704,508369,511693,513167,
        utr3 = self.genome.transcripts["NM_001083613"].utr3s[1] # TMEM219 gene
        self.assertEqual(Interval("chr16", "+", 513064, 513167, self.genome), utr3.interval)

        # partially coding exon (positive strand)
        tran = self.genome.transcripts["NM_001083613"]
        exon = tran.exons[4]
        self.assertEqual(exon.start, exon.cds.start)
        self.assertNotEqual(exon.start, exon.utr3.start)
        self.assertEqual(exon.end, exon.utr3.end)
        self.assertEqual(exon.cds.end, exon.utr3.start)

        # partially coding exon (negative strand)
        tran = self.genome.transcripts["NM_001109891"]
        exon = tran.exons[6]
        self.assertEqual(exon.start, exon.utr3.start)
        self.assertNotEqual(exon.end, exon.utr3.end)
        self.assertEqual(exon.end, exon.cds.end)
        self.assertEqual(exon.utr3.end, exon.cds.start)

        # non-coding exon
        tran = self.genome.transcripts["NM_001083613"]
        exon = tran.exons[5]
        self.assertIs(exon.cds, None)
        self.assertEqual(exon.interval, exon.utr3.interval)

        # UTR starting on exon boundary
        tran = self.genome.transcripts["NM_001310156"]
        self.assertEqual(1, len(tran.utr3s))
        self.assertEqual(tran.exons[-1].utr3.interval, tran.exons[-1].interval)
        self.assertIsNone(tran.exons[-1].cds)
        self.assertIsNone(tran.exons[-2].utr3)

########################################################

if __name__ == "__main__":
    unittest.main()
