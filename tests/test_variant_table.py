# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import unittest
import platform
import os
import gc

import numpy as np

from tempfile import mkstemp

from . import dumptext
from genome_kit import Interval
from genome_kit import VCFTable
from genome_kit import VCFVariant
from genome_kit import Variant
from . import MiniGenome


class TestVariant(unittest.TestCase):
    def test_attributes(self):

        # Create single interval one at a time.
        # Check initialization, getattr, repr, str, and del
        # Repeat several times to make sure allocation/refcounting works.
        while gc.collect():
            pass
        for _ in range(10):
            c = Variant("chr5", 2000, "AA", "T", "hg19")

            self.assertEqual(c.chrom, "chr5")
            self.assertEqual(c.chromosome, "chr5")
            self.assertEqual(c.strand, "+")
            self.assertTrue(c.is_positive_strand())
            self.assertEqual(c.start, 2000)
            self.assertEqual(c.end, 2002)
            self.assertEqual(c.end5, Interval("chr5", "+", 2000, 2000, "hg19"))
            self.assertEqual(c.end3, Interval("chr5", "+", 2002, 2002, "hg19"))
            self.assertEqual(c.refg, "hg19")
            self.assertEqual(c.sys, "hg19:chr5")
            self.assertIsInstance(c.__repr__(), str)
            self.assertIsInstance(c.__str__(), str)
            self.assertEqual(c.ref, "AA")
            self.assertEqual(c.alt, "T")

            # `anchor` and `anchor_offset` disabled on Variants.
            # This is the one way in which they CANNOT be treated as `Interval`
            with self.assertRaises(Exception):
                c.anchor
            with self.assertRaises(Exception):
                c.anchor_offset

            del c  # Make sure single dealloc works
            gc.collect()

    def test_refalt_smallstring(self):

        # Check that no bugs triggered when crossing ref/alt small string optimization threshold
        for n in range(10):
            for m in range(10):
                c = Variant("chr5", 2000, "A" * n, "T" * m, "hg19")
                self.assertEqual(c.ref, "A" * n)
                self.assertEqual(c.alt, "T" * m)

    def test_compare(self):

        v = Variant("chr1", 2000, "A", "T", "hg19")
        self.assertEqual(v, Variant("chr1", 2000, "A", "T", "hg19"))
        self.assertNotEqual(v, Variant("chr2", 2000, "A", "T", "hg19"))
        self.assertNotEqual(v, Variant("chr1", 2001, "A", "T", "hg19"))
        self.assertNotEqual(v, Variant("chr1", 2000, "C", "T", "hg19"))
        self.assertNotEqual(v, Variant("chr1", 2000, "A", "C", "hg19"))
        self.assertNotEqual(v, Variant("chr1", 2000, "A", "T", "hg38.p12"))

        w = Variant("chr1", 1000, "", "ATG", "hg19")
        self.assertEqual(w, Variant("chr1", 1000, "", "ATG", "hg19"))
        self.assertNotEqual(w, Variant("chr1", 1001, "", "ATG", "hg19"))
        self.assertIn(w, [Variant("chr1", 1000, "", "ATG", "hg19")])

    def test_from_variant_string(self):
        genome = MiniGenome("test_genome")
        self.assertEqual(Variant.from_string("chr1:13:A:T", genome), Variant("chr1", 12, "A", "T", genome))

        self.assertEqual(Variant.from_string("chr1:13:A:T", genome).as_variant_string(), "chr1:13:A:T")

        self.assertRaises(ValueError, Variant.from_string, "", genome)  # empty
        self.assertRaises(ValueError, Variant.from_string, "::", genome)  # empty and wrong num :
        self.assertRaises(ValueError, Variant.from_string, ":::", genome)  # empty and right num :
        self.assertRaises(ValueError, Variant.from_string, "chr:10:A:T", genome)  # chrom invalid
        self.assertRaises(ValueError, Variant.from_string, "chr1:10:A", genome)  # alt missing
        self.assertRaises(ValueError, Variant.from_string, "chr1:10:A:T:", genome)  # extra colon
        self.assertRaises(ValueError, Variant.from_string, "chr1:10:A:T:T", genome)  # extra field
        self.assertRaises(ValueError, Variant.from_string, "chr1::A:T", genome)  # start empty
        self.assertRaises(ValueError, Variant.from_string, "chr1:AT:A:T", genome)  # start not integer
        self.assertRaises(AttributeError, Variant.from_string, "chr1:10:TT:T", "hg38.p12")  # refg type invalid

    def test_sets(self):

        # If hash() or == for Variant is messed up, elements will be missing/repeated in sets.
        a = Variant("chr5", 2000, "A", "G", "hg19")
        b = Variant("chr5", 2000, "A", "G", "hg19")
        c = Variant("chr5", 2000, "C", "G", "hg19")
        d = Variant("chr5", 2000, "C", "G", "hg19")
        e = Variant("chr5", 2000, "A", "T", "hg19")
        f = Variant("chr5", 2000, "A", "T", "hg19")
        g = Variant("chr5", 2001, "A", "T", "hg19")
        h = Variant("chr5", 2001, "A", "T", "hg19")
        i = Variant("chr5", 2000, "A", "G", "hg38.p12")
        j = Variant("chr5", 2000, "A", "G", "hg38.p12")
        s1 = set([a, b, c, d, e, f, g, h, i, j])
        s2 = set([a, c, e, g, i])
        s3 = set([b, d, f, h, j])
        self.assertEqual(len(s1), 5)
        self.assertEqual(s1, s2)
        self.assertEqual(s1, s3)

    def test_hash(self):
        a = Variant("chr1", 155163643, "", "G", "hg19")
        b = Variant("chr1", 155163643, "", "G", "hg19")
        self.assertEqual(hash(a), hash(b))

    def test_serialize(self):

        variants = [
            Variant("chr1", 10, "", "", "hg19"),
            Variant("chr2", 11, "", "T", "hg38.p12"),
            Variant("chr1", 10, "A", "", "hg19"),
            Variant("chr1", 12, "A", "T", "hg19"),
            Variant("chrX", 12, "CAAAAAAAAAAAAAAAAAAAAAAAAAG", "T", "hg19"),
            Variant("chrY", 16, "A", "CTTTTTTTTTTTTTTTTTTTTTTTTTG", "hg19"),
        ]

        # Check that pickling / unpickling works
        import pickle
        stored = pickle.dumps(variants, pickle.HIGHEST_PROTOCOL)
        loaded = pickle.loads(stored)
        self.assertEqual(loaded, variants)

        # TODO test serialization of variants that come from VariantTable


#################################################################

vcf_header1 = """
    ##fileformat=VCFv4.2
    ##reference=file:///gk-local/hg19.fa
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##INFO=<ID=AD,Number=R,Type=Integer,Description="Total read depth for each allele">
    ##INFO=<ID=CIGAR,Number=A,Type=String,Description="Cigar string describing how to align an alternate allele to the reference allele">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
    ##INFO=<ID=POOLED,Number=.,Type=String,Description="Test multidimensional">
    ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
    ##INFO=<ID=VARLEN,Number=.,Type=String,Description="Test multidimensional">
    ##FORMAT=<ID=BADFLAG,Number=0,Type=Flag,Description="Not allowed by VCF specification">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
    ##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
    ##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allelic depths on the forward strand">
    ##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allelic depths on the reverse strand">
    ##FORMAT=<ID=EC,Number=A,Type=Integer,Description="Expected alternate allele counts">
    ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype quality">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">
    ##FORMAT=<ID=TESTSTR,Number=2,Type=String,Description="Test multidimensional">
    """

vcf_body1 = """
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SRR810007.splice.sorted.bam	SRR810031.splice.sorted.bam	SRR810051.splice.sorted.bam
    1	15015	.	G	C	222	.	AF=1;DP=228;OTHER_INFO_FIELDS=SKIPPED	GT:PL:DP:SP:ADF:ADR:AD	0/0:255,0,255:215:38:91,26:52,46:143,72	0/1:197,0,255:189:23:84,17:59,25:143,42	0/1:255,0,255:121:23:41,12:36,29:77,41
    1	15015	.	G	C	222	.	AF=1;DP=228;OTHER_INFO_FIELDS=SKIPPED	GT:PL:DP:SP:ADF:ADR:AD	0/1:255,0,255:215:38:91,26:52,46:143,72	0/1:197,0,255:189:23:84,17:59,25:143,42	0/1:255,0,255:121:23:41,12:36,29:77,41
    """


class TestVCFVariant(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

    def test_invalid_init(self):
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            VCFVariant('chr1', '+', 0, 0, MiniGenome())


def make_vcfbin(outfile, infile, reference_genome=MiniGenome('hg19'), validate=False, **kwargs):
    VCFTable.build_vcfbin(outfile, infile, reference_genome, validate=validate, **kwargs)
    return reference_genome


# If normalize=True, then from_vcf tries to dump a FASTA file of the genome DNA sequence
def from_vcf(vcfpath, reference_genome=MiniGenome('hg19'), validate=False, normalize=False, cache=False, **kwargs):
    return VCFTable.from_vcf(vcfpath, reference_genome, validate=validate, normalize=normalize, cache=cache, **kwargs)


class TestVCFTable(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        # Create a tmpfile that the test will continually clobber as it builds various tables
        fd, cls.tmpvcf = mkstemp(suffix=".vcf")
        os.close(fd)
        fd, cls.tmpbin = mkstemp(suffix=".vcf.bin")
        os.close(fd)
        cls.tmpgz = cls.tmpvcf + ".gz"
        cls.tmptbi = cls.tmpvcf + ".gz.tbi"

    @classmethod
    def tearDownClass(cls):  # Remove temporary files
        if os.path.exists(cls.tmpvcf): os.unlink(cls.tmpvcf)
        if os.path.exists(cls.tmpbin): os.unlink(cls.tmpbin)
        if os.path.exists(cls.tmpgz): os.unlink(cls.tmpgz)
        if os.path.exists(cls.tmptbi): os.unlink(cls.tmptbi)

    def make_vcfbin(self, reference_genome=MiniGenome('hg19'), validate=False, **kwargs):
        return make_vcfbin(self.tmpbin, self.tmpvcf, reference_genome, validate=validate, **kwargs)

    def from_vcf(self, reference_genome=MiniGenome('hg19'), validate=False, normalize=False, cache=False, **kwargs):
        return from_vcf(self.tmpvcf, reference_genome, validate=validate, normalize=normalize, cache=cache, vcfbinpath=self.tmpbin, **kwargs)

    def check_build_vcfbin(self, exception_type, msg, refg=MiniGenome("hg19"), **kwargs):
        with self.assertRaisesRegex(exception_type, msg):
            # TODO: VCFTable should allow mock Genomes
            # The VCFTable builder passes in the refg string from the python layer and
            # reconstitutes it as a C++ Genome object. Unfortunately that means you cannot create a
            # VCFTable from a mocked genome (eg. test genomes or GenBank ones).
            self.make_vcfbin(refg, **kwargs)

    def test_invalid_file(self):

        self.assertRaises(IOError, VCFTable, "non_existent_file_943412.vcf.bin")
        dumptext(self.tmpbin, "not_a_valid_vcfbin_file_header")
        self.assertRaises(IOError, VCFTable, self.tmpbin)

    def test_invalid_reference_genome(self):

        dumptext(self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            """)
        self.check_build_vcfbin(TypeError, "Genome", refg="h31")
        self.check_build_vcfbin(TypeError, "Genome", refg=37)

    def test_invalid_genotype(self):
        # GT requires the form .|. or ./. for missing
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT	..
            """)
        self.check_build_vcfbin(ValueError, "haploid.*diploid", fmt_ids=["GT"])

        # GT only accepts [0-9] or . for a value
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT	?/?
            """)
        self.check_build_vcfbin(ValueError, "GT component", fmt_ids=["GT"])

    def test_svtype(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	.
            1	123	.	G	C	123	.	SVTYPE=DEL
            1	123	.	G	C	123	.	SVTYPE=INS
            1	123	.	G	C	123	.	SVTYPE=DUP
            1	123	.	G	C	123	.	SVTYPE=INV
            1	123	.	G	C	123	.	SVTYPE=CNV
            1	123	.	G	C	123	.	SVTYPE=BND
            """)
        self.make_vcfbin(info_ids=["SVTYPE"])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].SVTYPE, VCFTable.SVTYPE_NA)
            self.assertEqual(vt[1].SVTYPE, VCFTable.SVTYPE_DEL)
            self.assertEqual(vt[2].SVTYPE, VCFTable.SVTYPE_INS)
            self.assertEqual(vt[3].SVTYPE, VCFTable.SVTYPE_DUP)
            self.assertEqual(vt[4].SVTYPE, VCFTable.SVTYPE_INV)
            self.assertEqual(vt[5].SVTYPE, VCFTable.SVTYPE_CNV)
            self.assertEqual(vt[6].SVTYPE, VCFTable.SVTYPE_BND)

    def test_id(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	.
            1	123	A	G	C	123	.	.
            1	123	B	G	C	123	.	.
            """)
        self.make_vcfbin(info_ids=["ID"])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].ID, "")
            self.assertEqual(vt[1].ID, "A")
            self.assertEqual(vt[2].ID, "B")

    def test_close(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1,1
            """)
        self.make_vcfbin(info_ids=["AF"], fmt_ids=["AD"])
        with VCFTable(self.tmpbin) as vt:
            variant = vt[0]
            refg = variant.refg

        interval = Interval("chr1", "+", 0, 0, refg)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_5p_aligned(interval)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_3p_aligned(interval)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_5p_within(interval)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_3p_within(interval)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_within(interval)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_overlapping(interval)
        with self.assertRaisesRegex(IOError, 'close'):
            vt.find_exact(interval)

        with self.assertRaisesRegex(IOError, 'close'):
            variant.chromosome

        with VCFTable(self.tmpbin) as vt:
            af = vt.info('AF')
            ad = vt.format('AD')

        np.testing.assert_equal(af, [1])
        np.testing.assert_equal(ad, [[[1, 1]]])
        af._mmap.close()
        ad._mmap.close()

    def test_value_count(self):
        # R == 2
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1,1
            """)
        self.make_vcfbin(fmt_ids=["AD"])
        # too few R
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1
            """)
        self.check_build_vcfbin(ValueError, "Expected 2 values, found 1", fmt_ids=["AD"])

        # too many R
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1,1,1
            """)
        self.check_build_vcfbin(ValueError, "Expected 2 values, found 3", fmt_ids=["AD"])

        # A == 1
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:EC	./.:1
            """)
        self.make_vcfbin(fmt_ids=["EC"])
        # too many A
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:EC	./.:1,1
            """)
        self.check_build_vcfbin(ValueError, "Expected 1 values, found 2", fmt_ids=["EC"])

        # G == 2
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:PL	./.:1,1
            """)
        self.make_vcfbin(fmt_ids=["PL"])
        # too few G
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:PL	./.:1
            """)
        self.check_build_vcfbin(ValueError, "Expected 2 values, found 1", fmt_ids=["PL"])
        # too many G
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:PL	./.:1,1,1
            """)
        self.check_build_vcfbin(ValueError, "Expected 2 values, found 3", fmt_ids=["PL"])

        # Number=2
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:HQ	./.:1,1
            """)
        self.make_vcfbin(fmt_ids=["HQ"])
        # too few 2
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:HQ	./.:1
            """)
        self.check_build_vcfbin(ValueError, "Expected 2 values, found 1", fmt_ids=["HQ"])
        # too many 2
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:HQ	./.:1,1,1
            """)
        self.check_build_vcfbin(ValueError, "Expected 2 values, found 3", fmt_ids=["HQ"])

        # Number=0
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	IMPRECISE
            """)
        self.make_vcfbin(info_ids=["IMPRECISE"])
        # too many 0
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	IMPRECISE=1
            """)
        self.check_build_vcfbin(ValueError, "Expected 0 values", info_ids=["IMPRECISE"])

    def test_variants_and_samples(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            1	123	.	G	C	123	.	AF=1;DP=123	GT:ADR:AD	0/0:123,123:123,1	0/1:123,123:132,2	1/1:123,123:123,3
            1	456	.	A	T	456	.	AF=1;DP=456	GT:ADR:AD	0/.:456,456:456,4	./1:456,456:465,5	./.:456,456:456,6
            """)
        self.make_vcfbin(fmt_ids=["GT", "AD"])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(len(vt), 2)
            self.assertEqual(vt.num_samples, 3)
            np.testing.assert_equal(vt.format("AD"), [[[123, 1], [132, 2], [123, 3]], [[456, 4], [465, 5], [456, 6]]])
            np.testing.assert_equal(vt.format("GT"), [[0, 1, 2], [3, 3, 3]])

            # Check invalid keys. Junk keys (AF, AL, FF) and omitted-from-build keys (DP)
            for k in ["ADR", "AF", "AL", "FF", "DP"]:
                self.assertRaises(KeyError, vt.format, k)

    def test_chrMT(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            MT	1	.	A	T	.	.	.
            chrM	2	.	A	T	.	.	.
            """)
        self.make_vcfbin(MiniGenome("hg38.p12"))
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(len(vt), 2)

    def test_exclude(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	1	.	A	T	.	.	.
            1	2	.	A	T	.	.	.
            3	3	.	A	T	.	.	.
            MT	4	.	A	T	.	.	.
            """)
        genome = MiniGenome("hg38.p12")
        self.make_vcfbin(genome,
            exclude=[Interval("chr1", "+", 1, 2, genome),
                       Interval("chrM", "+", 0, 1000, genome)])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].interval, Interval("chr1", "+", 0, 1, genome))
            self.assertEqual(vt[1].interval, Interval("chr3", "+", 2, 3, genome))

    def test_allow(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	1	.	A	T	.	.	.
            1	2	.	A	T	.	.	.
            3	3	.	A	T	.	.	.
            MT	4	.	A	T	.	.	.
            """)
        genome = MiniGenome("hg38.p12")
        self.make_vcfbin(genome,
            allow=[Interval("chr1", "+", 1, 2, genome),
                       Interval("chrM", "+", 0, 1000, genome)])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].interval, Interval("chr1", "+", 1, 2, genome))
            self.assertEqual(vt[1].interval, Interval("chrM", "+", 3, 4, genome))

    def test_variants_and_no_samples(self):

        # Check FORMAT but zero sample columns
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
            1	123	.	G	C	123	.	AF=1;DP=123	GT:ADR:AD
            1	456	.	A	T	456	.	AF=1;DP=456	GT:ADR:AD
            """)
        self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(len(vt), 2)
            self.assertEqual(vt.num_samples, 0)

        self.check_build_vcfbin(ValueError, "No FORMAT fields*", fmt_ids=["GT", "AD"])

        # Check no FORMAT column at all
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	AF=1;DP=123
            1	456	.	A	T	456	.	AF=1;DP=456
            """)
        self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(len(vt), 2)
            self.assertEqual(vt.num_samples, 0)

        self.check_build_vcfbin(ValueError, "No FORMAT fields", fmt_ids=["GT", "AD"])

    def test_samples_and_no_variants(self):

        dumptext(self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            """)
        self.make_vcfbin(fmt_ids=["GT", "AD"])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(len(vt), 0)
            self.assertEqual(vt.num_samples, 3)
            self.assertEqual(vt.format("AD").size, 0)
            self.assertEqual(vt.format("GT").size, 0)

    def test_no_samples_and_no_variants(self):

        dumptext(self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
            """)
        self.make_vcfbin(fmt_ids=["GT", "AD"])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(len(vt), 0)
            self.assertEqual(vt.num_samples, 0)
            self.assertEqual(vt.format("AD").size, 0)
            self.assertEqual(vt.format("GT").size, 0)

    def test_bad_header_line(self):

        # mandatory >= 8 cols
        dumptext(self.tmpvcf, vcf_header1, "#CHROM	POS	ID	REF	ALT	QUAL")
        self.check_build_vcfbin(ValueError, "Expected at least 8 tab")

        dumptext(self.tmpvcf, vcf_header1, "#CHROM	ALT	QUAL")
        self.check_build_vcfbin(ValueError, "Expected at least 8 tab")

    def test_parse_body_lines(self):

        # samples in header line must be consistent with subsequent rows
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            1	123	.	G	C	123	.	AF=1;DP=123	GT:ADR:AD	0/0:123,123:123,1	0/1:123,123:132,2
            """)
        self.check_build_vcfbin(ValueError, "Expects 3 samples but found 2", fmt_ids=['GT'])

        # Check specific REF/ALT combinations
        def check_refalt(ref, alt, err_str):
            dumptext(
                self.tmpvcf, vcf_header1, """
                #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
                1	123	.	%s	%s	123	.	AF=1;DP=123
                """ % (ref, alt))
            self.check_build_vcfbin(ValueError, err_str)

        check_refalt("", "A", "normalized")  # REF and ALT must be left aligned normalized
        check_refalt("G", "", "normalized")  # REF and ALT must be left aligned normalized
        check_refalt("-", "A", "Invalid REF")  # REF must be ACTGN
        check_refalt(".", "A", "Invalid REF")  # REF must be ACTGN
        check_refalt("G", "-", "Invalid ALT")  # ALT must be ACGTN.* or <SV>
        check_refalt("X", "A", "Invalid REF")  # REF must be "A", "C", "G", "T", or "N"
        check_refalt("G", "X", "Invalid ALT")  # ALT must be "A", "C", "G", "T", or "N"

    def test_normalize(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	122	.	CG	C	123	.	AF=1;DP=123
            1	122	.	C	CG	123	.	AF=1;DP=123
            """)
        self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].as_variant_string(), "chr1:123:G:")
            self.assertEqual(vt[1].as_variant_string(), "chr1:123::G")
            self.assertNotEqual(vt[0], vt[1])

    def test_missing(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3	samp4
            1	123	.	G	C	123	.	.	GT:AD:PL	./.	./.:.	./.:.:.	./.:.:1,1
        """)

        num_samples = 4

        # default zero fill
        self.make_vcfbin(
            info_ids=["IMPRECISE", "CIGAR", "AD"],
            fmt_ids=["GT", "PL"])
        with VCFTable(self.tmpbin) as vt:
            self.assertFalse(vt[0].IMPRECISE)
            self.assertEqual(vt[0].CIGAR, "")
            np.testing.assert_equal(vt.info("AD"), [[0, 0]])

            np.testing.assert_equal(vt.format("GT"), [num_samples * [VCFTable.GT_UNKNOWN]])
            np.testing.assert_equal(vt.format("PL"), [(num_samples - 1) * [2 * [0]] + [[1, 1]]])

        # user-specified fill value
        self.make_vcfbin(
            info_ids={"AD": 1},
            fmt_ids={
                "GT": None,
                "AD": None,
                "PL": -(2**31)
            })
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(vt.info("AD"), [[1, 1]])

            np.testing.assert_equal(vt.format("GT"), [num_samples * [VCFTable.GT_UNKNOWN]])
            np.testing.assert_equal(vt.format("AD"), [num_samples * [2 * [0]]])
            np.testing.assert_equal(vt.format("PL"), [(num_samples - 1) * [2 * [-(2**31)]] + [[1, 1]]])

    def test_undefined_field(self):
        # check if name matches but no header
        dumptext(
            self.tmpvcf, "", """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1.1;IMPRECISE;CIGAR=N	AD	1,1
        """)

        self.check_build_vcfbin(ValueError, "Missing value for INFO ID IMPRECISE", info_ids=["IMPRECISE"])
        self.check_build_vcfbin(ValueError, "header for ID CIGAR", info_ids=["CIGAR"])

        self.check_build_vcfbin(ValueError, "header for ID AF", info_ids=["AF"])
        self.check_build_vcfbin(ValueError, "Failed to parse", info_ids={"AF": np.int32(0)})
        self.make_vcfbin(info_ids={"AF": np.float32(0.0)})
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(vt.info("AF"), np.array([1.1], np.float32))

        self.check_build_vcfbin(ValueError, "header for ID AD", fmt_ids=["AD"])
        self.make_vcfbin(fmt_ids={"AD": np.int32(0)})
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(vt.format("AD"), [[[1, 1]]])
        self.make_vcfbin(fmt_ids={"AD": np.float32(0.0)})
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(vt.format("AD"), [[[1.0, 1.0]]])

        # check if name mismatches no header
        dumptext(self.tmpvcf, "", """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	DP=1
        """)
        self.check_build_vcfbin(ValueError, "header for ID IMPRECISE", info_ids=["IMPRECISE"])
        self.check_build_vcfbin(ValueError, "header for ID CIGAR", info_ids=["CIGAR"])

        self.check_build_vcfbin(ValueError, "header for ID AF", info_ids=["AF"])
        self.check_build_vcfbin(ValueError, "header for ID AF", info_ids={"AF": np.int32(0)})
        self.check_build_vcfbin(ValueError, "header for ID AF", info_ids={"AF": np.float32(0.0)})

        self.check_build_vcfbin(ValueError, "No FORMAT fields", fmt_ids=["AD"])
        self.check_build_vcfbin(ValueError, "No FORMAT fields", fmt_ids={"AD": np.int32(0)})
        self.check_build_vcfbin(ValueError, "No FORMAT fields", fmt_ids={"AD": np.float32(0.0)})

        # check if name mismatches with header
        dumptext(self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	DP=1
        """)

        self.make_vcfbin(info_ids=["AF", "AD", "IMPRECISE", "CIGAR"])
        with VCFTable(self.tmpbin) as vt:
            self.assertFalse(vt[0].IMPRECISE)
            self.assertEqual(vt[0].CIGAR, "")
            np.testing.assert_equal(vt.info("AD"), [[0, 0]])
            np.testing.assert_equal(vt.info("AF"), [0.0])

        self.check_build_vcfbin(ValueError, "AF is defined as a Float", info_ids={"AF": np.int32(0)})
        self.make_vcfbin(info_ids={"AF": np.float32(0.0)})
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(vt.info("AF"), [0.0])

        self.check_build_vcfbin(ValueError, "No FORMAT fields", fmt_ids=["AD"])
        self.check_build_vcfbin(ValueError, "No FORMAT fields", fmt_ids={"AD": np.int32(0)})
        self.check_build_vcfbin(ValueError, "AD defined as an Integer", fmt_ids={"AD": np.float32(0.0)})

        # check if format name mismatches with header
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1.1	GT	0/1
        """)
        with self.from_vcf(fmt_ids=["AD"]) as vt:
            np.testing.assert_equal(vt.format("AD"), [[[0, 0]]])
        with self.from_vcf(fmt_ids={"AD": np.int32(0)}) as vt:
            np.testing.assert_equal(vt.format("AD"), [[[0, 0]]])
        self.check_build_vcfbin(ValueError, "AD defined as an Integer", fmt_ids={"AD": np.float32(0.0)})

    def test_wrong_type(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	IMPRECISE;CIGAR=N;AF=1.1	AD	1,1
        """)

        self.check_build_vcfbin(ValueError, "AF is defined as a Float", info_ids={"AF": np.int8(0)})
        self.check_build_vcfbin(ValueError, "AF is defined as a Float", info_ids={"AF": np.int16(0)})
        self.check_build_vcfbin(ValueError, "AF is defined as a Float", info_ids={"AF": np.int32(0)})

        self.check_build_vcfbin(ValueError, "AD defined as an Integer", fmt_ids={"AD": np.float16(0.0)})
        self.check_build_vcfbin(ValueError, "AD defined as an Integer", fmt_ids={"AD": np.float32(0.0)})

        self.check_build_vcfbin(ValueError, "CIGAR is defined as a String", info_ids={"CIGAR": 1})

        self.check_build_vcfbin(ValueError, "IMPRECISE is defined as a Flag", info_ids={"IMPRECISE": 1})

    def test_float(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	C	123	.	AF=1.1
            1	123	.	G	C	123	.	AF=3.4028235E+38
            1	123	.	G	C	123	.	AF=-3.4028235e+38
            1	123	.	G	C	123	.	AF=+.1E-10
            1	123	.	G	C	123	.	AF=NaN
            1	123	.	G	C	123	.	AF=Inf
            1	123	.	G	C	123	.	AF=+Inf
            1	123	.	G	C	123	.	AF=-Inf
        """)
        self.make_vcfbin(info_ids=["AF"])
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(
                vt.info("AF"),
                np.array([
                    1.1, 3.4028235E+38, -3.4028235e+38, +.1E-10,
                    float('NaN'),
                    float('Inf'),
                    float('+Inf'),
                    float('-Inf')
                ], np.float32))

        self.check_build_vcfbin(ValueError, "Overflow", info_ids={"AF": np.float16(0.0)})

        dumptext(
            self.tmpvcf, vcf_header1, """
                #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
                1	123	.	G	C	123	.	AF=1.1
                1	123	.	G	C	123	.	AF=65504.0E0
                1	123	.	G	C	123	.	AF=-65504.0E0
                1	123	.	G	C	123	.	AF=+.1E-5
                1	123	.	G	C	123	.	AF=NaN
                1	123	.	G	C	123	.	AF=Inf
                1	123	.	G	C	123	.	AF=+Inf
                1	123	.	G	C	123	.	AF=-Inf
            """)
        self.make_vcfbin(info_ids={"AF": np.float16(0.0)})
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_almost_equal(
                vt.info("AF"),
                np.array([1.1, 65504.0E0, -65504.0E0, +.1E-5,
                          float('NaN'),
                          float('Inf'),
                          float('+Inf'),
                          float('-Inf')], np.float16), 2)

            self.assertIsInstance(vt[0].AF, float)
            np.testing.assert_almost_equal(vt[0].AF, 1.1, 2)
            np.testing.assert_almost_equal(vt[1].AF, 65504.0E0)
            np.testing.assert_almost_equal(vt[2].AF, -65504.0E0)
            np.testing.assert_almost_equal(vt[3].AF, +.1E-5)
            np.testing.assert_equal(vt[4].AF, float('NaN'))
            np.testing.assert_equal(vt[5].AF, float('Inf'))
            np.testing.assert_equal(vt[6].AF, float('+Inf'))
            np.testing.assert_equal(vt[7].AF, float('-Inf'))

    def test_str(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2
            1	123	.	G	C	123	.	CIGAR=1M;VARLEN=1M,%3A%3B%3D%25%2C%0D%0A%09,.	TESTSTR	1M,%3A%3B%3D%25%2C%0D%0A%09	.
            1	123	.	G	C	123	.	CIGAR=2M;VARLEN=0,1,1.1,NaN,Inf	TESTSTR	0,1	1.1,NaN
            1	123	.	G	C	123	.	CIGAR=3M;VARLEN=%3A%3B%3D%25%2C%0D%0A%09	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=4M;VARLEN=.	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=5M;VARLEN=.	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=6M;VARLEN=.	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=7M;VARLEN=.	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=8M;VARLEN=.	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=9M;VARLEN=SecondLast	TESTSTR	.	.
            1	123	.	G	C	123	.	CIGAR=10M;VARLEN=Last	TESTSTR	.	.
        """)

        # str
        self.make_vcfbin(info_ids=["CIGAR", "VARLEN"])
        with VCFTable(self.tmpbin) as vt:
            with self.assertRaises(TypeError):
                vt.info("CIGAR")

            for i in range(0, len(vt)):
                self.assertEqual(vt[i].CIGAR, "{}M".format(i + 1))

            # no split no decode
            #self.assertListEqual(vt[0].VARLEN, ["1M", ":;=%,\r\n\t", ""])
            #self.assertListEqual(vt[1].VARLEN, ["0", "1", "1.1", "NaN", "Inf"])
            #self.assertEqual(vt[0].VARLEN, ":;=%,\r\n\t")
            self.assertEqual(vt[0].VARLEN, "1M,%3A%3B%3D%25%2C%0D%0A%09,.")
            self.assertEqual(vt[1].VARLEN, "0,1,1.1,NaN,Inf")
            self.assertEqual(vt[2].VARLEN, "%3A%3B%3D%25%2C%0D%0A%09")
            for i in range(3, len(vt) - 2):
                self.assertEqual(vt[i].VARLEN, "")
            self.assertEqual(vt[-2].VARLEN, "SecondLast")
            self.assertEqual(vt[-1].VARLEN, "Last")

            with self.assertRaisesRegex(TypeError, "VCFVariant attribute"):
                vt.info("CIGAR")

        # FORMAT unsupported
        self.check_build_vcfbin(ValueError, "String type.*unsupported", fmt_ids=["TESTSTR"])

    def test_flag(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	IMPRECISE	BADFLAG	.
            1	123	.	G	C	123	.	.	.	.
            1	123	.	G	C	123	.	AD=0	.	.
            """)

        self.make_vcfbin(info_ids=['IMPRECISE'])
        with VCFTable(self.tmpbin) as vt:
            self.assertTrue(vt[0].IMPRECISE)
            self.assertFalse(vt[1].IMPRECISE)
            self.assertFalse(vt[2].IMPRECISE)

        self.check_build_vcfbin(ValueError, "FLAG types.*unsupported", fmt_ids=["BADFLAG"])

    def test_end_svtypes(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	.	G	<DEL>	123	.	SVTYPE=DEL;END=130
            1	123	.	G	<INS>	123	.	SVTYPE=INS;END=130
            1	123	.	G	<DUP>	123	.	SVTYPE=DUP;END=130
            1	123	.	G	<INV>	123	.	SVTYPE=INV;END=130
            1	123	.	G	<CNV>	123	.	SVTYPE=CNV;END=130
            1	123	.	G	<DUP:TANDEM>	123	.	SVTYPE=DUP;END=130
            1	123	.	G	<DEL:ME>	123	.	SVTYPE=DEL;END=130
            1	123	.	G	<INS:ME>	123	.	SVTYPE=INS;END=130
            """)

        genome = self.make_vcfbin(info_ids=['SVTYPE'])
        with VCFTable(self.tmpbin) as vt:
            np.testing.assert_equal(
                vt.info('SVTYPE'), [
                    VCFTable.SVTYPE_DEL,
                    VCFTable.SVTYPE_INS,
                    VCFTable.SVTYPE_DUP,
                    VCFTable.SVTYPE_INV,
                    VCFTable.SVTYPE_CNV,
                    VCFTable.SVTYPE_DUP,
                    VCFTable.SVTYPE_DEL,
                    VCFTable.SVTYPE_INS,
                ])
            self.assertEqual([x.interval for x in vt], [Interval('chr1', '+', 123, 130, genome)] * 8)
            self.assertEqual([x.alt for x in vt],
                             ["<DEL>", "<INS>", "<DUP>", "<INV>", "<CNV>", "<DUP:TANDEM>", "<DEL:ME>", "<INS:ME>"])

            with self.assertRaises(ValueError):
                genome.variant_dna(vt[0].interval, vt[0])

    def test_svtype_bnd(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	right_join_after	G	G[2:123[	123	.	SVTYPE=BND
            1	123	left_join_after	G	G]2:123]	123	.	SVTYPE=BND
            1	123	left_join_before	G	]2:123]G	123	.	SVTYPE=BND
            1	123	right_join_before	G	[2:123[G	123	.	SVTYPE=BND
            1	123	insert_before	G	]2:123]NNNNG	123	.	SVTYPE=BND
            1	123	insert_after	G	GNNNN[2:123[	123	.	SVTYPE=BND
            """)
        genome = self.make_vcfbin(info_ids=['SVTYPE'])
        after_interval = Interval('chr1', '+', 123, 123, genome)
        before_interval = after_interval.shift(-1)
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].interval, after_interval)
            self.assertEqual(vt[1].interval, after_interval)
            self.assertEqual(vt[2].interval, before_interval)
            self.assertEqual(vt[3].interval, before_interval)
            self.assertEqual(vt[4].interval, before_interval)
            self.assertEqual(vt[5].interval, after_interval)

            self.assertEqual([x.ref for x in vt], [""] * 6)
            self.assertEqual(vt[0].alt, "[2:123[")
            self.assertEqual(vt[1].alt, "]2:123]")
            self.assertEqual(vt[2].alt, "]2:123]")
            self.assertEqual(vt[3].alt, "[2:123[")
            self.assertEqual(vt[4].alt, "]2:123]NNNN")
            self.assertEqual(vt[5].alt, "NNNN[2:123[")

    def test_svtype_bnd_multi_mate(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	123	multi	G	G[2:123[,G[3:123[	123	.	SVTYPE=BND
            """)
        self.check_build_vcfbin(ValueError, "multiallelic")

    @unittest.expectedFailure
    def test_svtype_bnd_telomere(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	0	telomere	N	.[2:124[	123	.	SVTYPE=BND
            1	1	multi	G	]2:123]G	123	.	SVTYPE=BND
            """)
        genome = self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt[0].interval, Interval('chr1', '+', -1, -2, genome))
            self.assertEqual(vt[1].interval, Interval('chr1', '+', 0, -1, genome))

    def test_readonly_error(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD:PL	./.	./.:.
        """)

        # Check tables are read-only
        self.make_vcfbin(fmt_ids=["GT", "PL"])
        with VCFTable(self.tmpbin) as vt:
            with self.assertRaisesRegex(ValueError, "read-only"):
                vt.format("GT")[:] = 0
            with self.assertRaisesRegex(ValueError, "read-only"):
                vt.format("PL")[:] = 0

    def test_format_dtypes(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD:PL	./.	./.:4,5
        """)

        # Check default dtypes
        self.make_vcfbin(fmt_ids=["GT", "AD"])
        with VCFTable(self.tmpbin) as vt:
            self.assertEqual(vt.format("GT").dtype, np.int8)
            self.assertEqual(vt.format("AD").dtype, np.int32)

        # Check custom dtypes
        fmt_dict = {"GT": VCFTable.GT_UNKNOWN, "AD": np.int16(-1)}
        self.make_vcfbin(fmt_ids=fmt_dict)
        with VCFTable(self.tmpbin) as vt:
            gt = vt.format("GT")
            ad = vt.format("AD")
            self.assertEqual(gt.dtype, np.int8)
            self.assertEqual(ad.dtype, np.int16)
            np.testing.assert_equal(gt[0, 0], VCFTable.GT_UNKNOWN)
            np.testing.assert_equal(ad[0, 0], [-1, -1])

    def test_info(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	122	.	AGC	A	123	.	AD=50,73;DP=123	AD	2,3
            1	456	.	A	TAC	456	.	AD=.,256;.	AD	4,5
            1	456	.	A	TAC	456	.	.	AD	6,7
            """)
        self.make_vcfbin(info_ids=["AD", "DP"], fmt_ids=["AD"])
        with VCFTable(self.tmpbin) as vt:
            ad = vt.info("AD")
            dp = vt.info("DP")
            ad_format = vt.format("AD")
            np.testing.assert_equal(ad, [[50, 73], [0, 256], [0, 0]])
            np.testing.assert_equal(dp, [123, 0, 0])
            np.testing.assert_equal(ad_format, [[[2, 3]], [[4, 5]], [[6, 7]]])

            with self.assertRaisesRegex(ValueError, "scalar values"):
                vt[0].AD

    def test_table(self):

        # Check basic table with 3 variants, two unique
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            1	122	.	AGC	A	123	.	AF=1;DP=123	GT:ADR:AD	0/0:123,123:123,1	0/1:123,123:132,2	1/1:123,123:123,3
            1	456	.	A	TAC	456	.	AF=1;DP=456	GT:ADR:AD	0/0:456,456:456,4	0/1:456,456:465,5	./.:456,456:456,6
            1	456	.	A	TAC	456	.	AF=1;DP=456	GT:ADR:AD	0/0:456,456:456,4	0/1:456,456:465,5	./.:456,456:456,6
            """)
        genome = self.make_vcfbin(fmt_ids=["AD"])
        with VCFTable(self.tmpbin) as vt:
            # Check length
            self.assertEqual(len(vt), 3)

            # Check variant 0
            v0 = vt[0]
            self.assertIsInstance(v0, VCFVariant)
            self.assertIsInstance(v0, Variant)
            self.assertEqual(v0.start, 122)
            self.assertEqual(v0.strand, "+")
            self.assertEqual(v0.ref, 'GC')
            self.assertEqual(v0.alt, '')
            self.assertEqual(len(v0), 2)
            self.assertEqual(v0.refg, genome.refg)
            self.assertEqual(v0.as_variant_string(), 'chr1:123:GC:')

            # Check variant 1
            v1 = vt[1]
            self.assertEqual(v1.start, 455)
            self.assertEqual(v1.strand, "+")
            self.assertEqual(v1.ref, 'A')
            self.assertEqual(v1.alt, 'TAC')
            self.assertEqual(len(v1), 1)
            self.assertEqual(v1.refg, genome.refg)
            self.assertEqual(v1.as_variant_string(), 'chr1:456:A:TAC')

            # Check variant 2
            v2 = vt[2]
            self.assertEqual(v1, v2)
            self.assertEqual(v2.as_variant_string(), v1.as_variant_string())

            # Check equality
            self.assertEqual(v0, v0)
            self.assertEqual(v1, v1)
            self.assertNotEqual(v0, v1)
            self.assertNotEqual(v0, v2)

            # Check range queries
            self.assertEqual(vt.find_5p_aligned(v0.end5), [v0])
            self.assertEqual(vt.find_5p_aligned(v1.end5), [v1, v2])
            self.assertEqual(vt.find_3p_aligned(v0.end3), [v0])
            self.assertEqual(vt.find_3p_aligned(v1.end3), [v1, v2])
            self.assertEqual(vt.find_exact(v0), [v0])
            self.assertEqual(vt.find_exact(v1), [v1, v2])
            self.assertEqual(vt.find_within(v0.expand(0, 500)), [v0, v1, v2])
            self.assertEqual(vt.find_overlapping(v0.shift(1)), [v0])
            self.assertEqual(vt.find_overlapping(v0.shift(3)), [])

            # Ensure table is unstranded
            self.assertEqual(vt.find_exact(v1.as_negative_strand()), [v1, v2])
            self.assertEqual(vt.find_within(v0.expand(0, 500).as_negative_strand()), [v0, v1, v2])
            self.assertEqual(vt.find_overlapping(v0.shift(1).as_negative_strand()), [v0])
            with self.assertRaises(ValueError):
                vt.find_5p_aligned(v0.as_negative_strand())
            with self.assertRaises(ValueError):
                vt.find_3p_aligned(v0.as_negative_strand())
            with self.assertRaises(ValueError):
                vt.find_5p_within(v0.as_negative_strand())
            with self.assertRaises(ValueError):
                vt.find_3p_within(v0.as_negative_strand())

            # Check where() lookup by mask
            ad = vt.format("AD")
            self.assertEqual(vt.where(np.any(ad.sum(axis=2) == 126, axis=1)), [v0])
            self.assertEqual(vt.where(np.any(ad.sum(axis=2) == 470, axis=1)), [v1, v2])

    def test_table_variant_equality(self):

        # Check VCFVariant equality/inequality
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            1	20	.	A	T	1	.	AF=1;DP=123
            1	20	.	A	T	1	.	AF=1;DP=123
            1	20	.	A	C	1	.	AF=1;DP=456
            1	20	.	G	T	1	.	AF=1;DP=456
            """)
        self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            a, b, c, d = vt
            self.assertEqual(a, b)
            self.assertNotEqual(a, c)
            self.assertNotEqual(a, d)
            self.assertNotEqual(b, c)
            self.assertNotEqual(b, d)
            self.assertNotEqual(c, d)

    # TODO: replace CI with a conda-forge specific env var
    @unittest.skipIf("CI" in os.environ, "bioconda dependencies not supported on conda-forge")
    def test_from_vcf(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
            2	20	.	A	T	1	.	AF=1;DP=12
            2	20	.	A	T	1	.	AF=1;DP=34
            2	20	.	A	C	1	.	AF=1;DP=56
            2	21	.	AG	ATC	1	.	AF=1;DP=78
            """)

        is_windows = (platform.system() == "Windows")
        if is_windows:
            tmpvcf = self.tmpvcf
        else:
            self.assertFalse(os.system("bgzip " + self.tmpvcf))
            self.assertFalse(os.system("tabix -p vcf -f " + self.tmpgz))
            tmpvcf = self.tmpgz


        os.unlink(self.tmpbin)  # Ensure binary file is gone
        self.assertFalse(os.path.exists(self.tmpbin))

        # No normalize, no cache
        with from_vcf(tmpvcf, info_ids=["AF", "DP"], vcfbinpath=self.tmpbin) as vt:
            self.assertEqual(len(vt), 4)
            self.assertEqual(vt[2].DP, 56)

        self.assertTrue(os.path.exists(self.tmpbin))

        # Cache
        with from_vcf(tmpvcf, info_ids=["AF", "DP"], vcfbinpath=self.tmpbin) as vt:
            self.assertEqual(len(vt), 4)
            self.assertEqual(vt[2].DP, 56)

        # Normalize
        if not is_windows:
            with from_vcf(tmpvcf, info_ids=["AF", "DP"], vcfbinpath=self.tmpbin) as vt:
                self.assertEqual(len(vt), 4)
                self.assertEqual(vt[3].ref, 'G')
                self.assertEqual(vt[3].alt, 'TC')
                self.assertEqual(vt[3].DP, 78)

    def test_info_ids(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1,1
            """)
        for ids in [[], ["AF"], ["DP"], ["AF", "DP"]]:
            with self.from_vcf(info_ids=ids) as vt:
                self.assertEqual(sorted(vt.info_ids), ids)

    def test_format_ids(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1,1
            """)
        for ids in [[], ["GT"], ["AD"], ["AD", "GT"]]:
            with self.from_vcf(fmt_ids=ids) as vt:
                self.assertEqual(sorted(vt.format_ids), ids)

    def test_sample_names(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            1	123	.	G	C	123	.	AF=1;DP=123	GT:AD	./.:1,1	./.:1,1	./.:1,1
            """)
        with self.from_vcf() as vt:
            self.assertEqual(vt.sample_names, [])

        with self.from_vcf(fmt_ids=["GT"]) as vt:
            self.assertEqual(vt.sample_names, ["samp1", "samp2", "samp3"])

    def test_ancenstral_variants(self):
        #  note ancestral variants always pass validation due to how GK trims the ref to empty
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            1	123	.	G	.	123	.	AF=1;DP=123	GT:AD	./.:1,1	./.:1,1	./.:1,1
            """)
        with self.assertRaisesRegex(ValueError, "Ancestral"):
            self.make_vcfbin()
        with self.from_vcf(ancestral="warn") as vt:
            self.assertEqual(str(vt[0]), "chr1:124::")
        with self.from_vcf(ancestral="exclude") as vt:
            self.assertEqual(len(vt), 0)

        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1	samp2	samp3
            1	123	.	G	G	123	.	AF=1;DP=123	GT:AD	./.:1,1	./.:1,1	./.:1,1
            """ )
        with self.assertRaisesRegex(ValueError, "alternative allele"):
            self.make_vcfbin()

    def test_sequence_variations_no_hit(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	AGC	A	.	.	AF=0.5	GT	0/1
            """)

        expected = []

        genome = self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            vcf_variants = [(x.chromosome, x.start, x.ref, x.alt)
                            for x in vt.sequence_variations(Interval("chr1", "+", 95, 95, genome, 95))]
        self.assertEqual(vcf_variants, expected)

    def test_sequence_variations_sorted(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	102	.	A	AT	.	.	AF=0.5	GT	0/1
            1	100	.	AGC	A	.	.	AF=0.5	GT	0/1
            1	104	.	AGC	AT	.	.	AF=0.5	GT	0/1
            """)

        expected = [("chr1", 100, "GC", ""), ("chr1", 102, "", "T"), ("chr1", 104, "GC", "T")]

        genome = self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            vcf_variants = [(x.chromosome, x.start, x.ref, x.alt)
                            for x in vt.sequence_variations(Interval("chr1", "+", 100, 105, genome))]
        self.assertEqual(vcf_variants, expected)

    def test_sequence_variations_insert_at_pos5(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	AGC	A	.	.	AF=0.5	GT	0/1
            1	102	.	A	AT	.	.	AF=0.5	GT	0/1
            1	104	.	AGC	AT	.	.	AF=0.5	GT	0/1
            """)

        expected = [("chr1", 102, "", "T"), ("chr1", 104, "GC", "T")]

        genome = self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            vcf_variants = [(x.chromosome, x.start, x.ref, x.alt)
                            for x in vt.sequence_variations(Interval("chr1", "+", 102, 105, genome))]
        self.assertEqual(vcf_variants, expected)

    def test_sequence_variations_only_insertions(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	A	AGC	.	.	AF=0.5	GT	0/1
            1	101	.	A	T	.	.	AF=0.5	GT	0/1
            """)

        expected = [("chr1", 100, "", "GC")]

        genome = self.make_vcfbin()
        with VCFTable(self.tmpbin) as vt:
            vcf_variants = [(x.chromosome, x.start, x.ref, x.alt)
                            for x in vt.sequence_variations(Interval("chr1", "+", 100, 100, genome))]
        self.assertEqual(vcf_variants, expected)

    def test_sequence_variations_anchor_right(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	102	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	104	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	106	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	108	.	AGC	A	.	.	AF=1.0	GT	1/1
            """)

        variants = [("chr1", 100, "GC", ""), ("chr1", 102, "GC", ""), ("chr1", 104, "GC", ""), ("chr1", 106, "GC", ""),
                    ("chr1", 108, "GC", "")]

        genome = self.make_vcfbin(info_ids=["AF"], fmt_ids=["GT"])
        with VCFTable(self.tmpbin) as vt:
            vcf_variants = [(x.chromosome, x.start, x.ref, x.alt)
                            for x in vt.sequence_variations(Interval("chr1", "+", 104, 106, genome, anchor='3p'))]
        self.assertEqual(vcf_variants, variants[:3])

    def test_sequence_variations_anchor_left(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	102	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	104	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	106	.	AGC	A	.	.	AF=1.0	GT	1/1
            1	108	.	AGC	A	.	.	AF=1.0	GT	1/1
            """)

        variants = [("chr1", 100, "GC", ""), ("chr1", 102, "GC", ""), ("chr1", 104, "GC", ""), ("chr1", 106, "GC", ""),
                    ("chr1", 108, "GC", "")]

        genome = self.make_vcfbin(info_ids=["AF"], fmt_ids=["GT"])
        with VCFTable(self.tmpbin) as vt:
            vcf_variants = [(x.chromosome, x.start, x.ref, x.alt)
                            for x in vt.sequence_variations(Interval("chr1", "+", 104, 106, genome, anchor='5p'))]
        self.assertEqual(vcf_variants, variants[2:])

    def test_masked(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	AGC	A	.	.	AF=0.1	GT	.
            1	102	.	AGC	A	.	.	AF=0.2	GT	0/0
            1	104	.	AGC	A	.	.	AF=0.3	GT	1/0
            1	106	.	AGC	A	.	.	AF=0.4	GT	1/1
            1	108	.	AGC	A	.	.	AF=0.5	GT	.
            """)


        genome = self.make_vcfbin(info_ids=["AF"], fmt_ids=["GT"])
        with VCFTable(self.tmpbin) as vt:
            with vt.masked([vt.index_of(x)
                            for x in vt.find_overlapping(genome.interval('chr1', '+', 102, 107))]) as index_masked:
                np.testing.assert_equal([x.start for x in index_masked], [102, 104, 106])
                np.testing.assert_equal(index_masked.info('AF'), np.array([0.2, 0.3, 0.4], dtype=np.float32))
                np.testing.assert_equal(
                    index_masked.format('GT'),
                    [[VCFTable.GT_HOMOZYGOUS_REF], [VCFTable.GT_HETEROZYGOUS_UNPHASED], [VCFTable.GT_HOMOZYGOUS_ALT]])
                with index_masked.masked(np.array([True, False, True])) as bitmasked:
                    np.testing.assert_equal([x.start for x in bitmasked], [102, 106])
                    np.testing.assert_equal(bitmasked.info('AF'), np.array([0.2, 0.4], dtype=np.float32))
                    np.testing.assert_equal(bitmasked.format('GT'),
                                            [[VCFTable.GT_HOMOZYGOUS_REF], [VCFTable.GT_HOMOZYGOUS_ALT]])

            # test if we do bitmask then index to see if we get the same results
            with vt.masked(np.array([True, True, False, True, True])) as bitmasked:
                np.testing.assert_equal([x.start for x in bitmasked], [100, 102, 106, 108])
                np.testing.assert_equal(bitmasked.info('AF'), np.array([0.1, 0.2, 0.4, 0.5], dtype=np.float32))
                np.testing.assert_equal(bitmasked.format('GT'), [[VCFTable.GT_UNKNOWN], [VCFTable.GT_HOMOZYGOUS_REF],
                                                                 [VCFTable.GT_HOMOZYGOUS_ALT], [VCFTable.GT_UNKNOWN]])
                with bitmasked.masked(
                    [bitmasked.index_of(x)
                     for x in bitmasked.find_overlapping(genome.interval('chr1', '+', 102, 107))]) as index_masked:
                    np.testing.assert_equal([x.start for x in index_masked], [102, 106])
                    np.testing.assert_equal(index_masked.info('AF'), np.array([0.2, 0.4], dtype=np.float32))
                    np.testing.assert_equal(index_masked.format('GT'),
                                            [[VCFTable.GT_HOMOZYGOUS_REF], [VCFTable.GT_HOMOZYGOUS_ALT]])

    def test_phased(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	100	.	AGC	A	.	.	AF=0.1	GT	.
            1	102	.	AGC	A	.	.	AF=0.2	GT	0/0
            1	104	.	AGC	A	.	.	AF=0.3	GT	1/0
            1	106	.	AGC	A	.	.	AF=0.3	GT	0/1
            1	108	.	AGC	A	.	.	AF=0.4	GT	1/1
            1	110	.	AGC	A	.	.	AF=0.5	GT	0|0
            1	112	.	AGC	A	.	.	AF=0.3	GT	1|0
            1	114	.	AGC	A	.	.	AF=0.3	GT	0|1
            1	116	.	AGC	A	.	.	AF=0.4	GT	1|1
            """)

        genome = MiniGenome("hg19")

        VCFTable.build_vcfbin(self.tmpbin, self.tmpvcf, genome, info_ids=["AF"], fmt_ids=["GT"], validate=False)
        with VCFTable(self.tmpbin) as vt:
            variants = vt.find_within(genome.interval('chr1', '+', 100, 120))
            self.assertEqual(VCFTable.GT_UNKNOWN, variants[0].GT)
            self.assertEqual(VCFTable.GT_HOMOZYGOUS_REF, variants[1].GT)
            self.assertEqual(VCFTable.GT_HETEROZYGOUS_UNPHASED, variants[2].GT)
            self.assertEqual(VCFTable.GT_HETEROZYGOUS_UNPHASED, variants[3].GT)
            self.assertEqual(VCFTable.GT_HOMOZYGOUS_ALT, variants[4].GT)
            self.assertEqual(VCFTable.GT_HOMOZYGOUS_REF, variants[5].GT)
            self.assertEqual(VCFTable.GT_HETEROZYGOUS_PHASED_1_0, variants[6].GT)
            self.assertEqual(VCFTable.GT_HETEROZYGOUS_PHASED_0_1, variants[7].GT)
            self.assertEqual(VCFTable.GT_HOMOZYGOUS_ALT, variants[8].GT)


    def test_phased_unknown(self):
        dumptext(
            self.tmpvcf, vcf_header1, """
            #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp1
            1	102	.	AGC	A	.	.	AF=0.2	GT	0/.
            1	104	.	AGC	A	.	.	AF=0.3	GT	1/.
            1	106	.	AGC	A	.	.	AF=0.4	GT	./0
            1	108	.	AGC	A	.	.	AF=0.5	GT	./1
            1	112	.	AGC	A	.	.	AF=0.7	GT	0|.
            1	114	.	AGC	A	.	.	AF=0.8	GT	1|.
            1	116	.	AGC	A	.	.	AF=0.9	GT	.|0
            1	118	.	AGC	A	.	.	AF=0.1	GT	.|1
            1	110	.	AGC	A	.	.	AF=0.6	GT	./.
            1	120	.	AGC	A	.	.	AF=0.2	GT	.|.
            """)

        genome = MiniGenome("hg19")

        VCFTable.build_vcfbin(self.tmpbin, self.tmpvcf, genome, info_ids=["AF"], fmt_ids=["GT"], validate=False)
        with VCFTable(self.tmpbin) as vt:
            variants = vt.find_within(genome.interval('chr1', '+', 100, 122))
            for i in range(10):
                self.assertEqual(VCFTable.GT_UNKNOWN, variants[i].GT)

        VCFTable.build_vcfbin(self.tmpbin, self.tmpvcf, genome, info_ids=["AF"], fmt_ids={"GT": VCFTable.GT_HETEROZYGOUS_UNPHASED}, validate=False)
        with VCFTable(self.tmpbin) as vt:
            variants = vt.find_within(genome.interval('chr1', '+', 100, 122))
            for i in range(4): # unphased, default value used
                self.assertEqual(VCFTable.GT_HETEROZYGOUS_UNPHASED, variants[i].GT)
            self.assertEqual(VCFTable.GT_UNKNOWN, variants[4].GT)
            for i in range(4, 8): # phased, default value ignored
                self.assertEqual(VCFTable.GT_UNKNOWN, variants[i].GT)
            for i in range(8, 10): # ./. and .|.
                self.assertEqual(VCFTable.GT_UNKNOWN, variants[i].GT)

if __name__ == "__main__":
    unittest.main()
