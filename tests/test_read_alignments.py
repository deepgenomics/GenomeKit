# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import unittest
import tempfile
import shutil
import os
import random
import time
import re
import numpy as np
from random import randint
from tempfile import mkstemp
from genome_kit import Alignment
from genome_kit import AlignmentMatch
from genome_kit import Interval
from genome_kit import Genome
from genome_kit import ReadAlignments
from genome_kit import JReadAlignments
from genome_kit import ReadDistributions
from genome_kit import Junction
from genome_kit import JunctionReadAlignmentsTable
from genome_kit import JunctionReadAlignments
from genome_kit import JunctionReadAlignment
from genome_kit import JunctionReadDistribution
from genome_kit import VariantTable
from genome_kit import Variant
from . import MiniGenome
from . import check_pythonic_indexing
from . import dumptext

TEST_DATA_FOLDER = os.path.join('tests', 'data', 'read_alignments')


def time_func(f):
    def time_function(*args, **kwargs):
        t0 = time.time()
        result = f(*args, **kwargs)
        t1 = time.time()
        print("Time Took %s: %s seconds" % (f.func_name, str(t1 - t0)))
        return result

    return time_function


class TempDir(object):
    def __init__(self):
        self.path = tempfile.mkdtemp()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if len(self.path) > 8:  # Try to avoid rm -r "/"
            shutil.rmtree(self.path)

    def __call__(self, relative_path):
        return os.path.join(self.path, relative_path)


def sample_junctions(num_junctions, seed=None):
    """pick out some junctions from annotations
    """

    if seed is not None:
        random.seed(seed)

    genome = MiniGenome()

    # the returned intervals should be unique and have non-zero length.
    junctions = list(set([intron.interval for intron in genome.introns if len(intron) > 0]))
    random.shuffle(junctions)
    return junctions[:num_junctions]


class SamFile(object):
    def __init__(self, filename, genome_version):
        self.filename = filename
        self.genome_version = genome_version
        self.samfile = None

    def __enter__(self):
        self.samfile = open(self.filename, 'w')
        self._write_header()
        return self

    def __exit__(self, typ, value, traceback):
        if self.samfile is not None:
            self.samfile.close()
            self.samfile = None

    def _write_header(self):
        header = '@PG	CL:"i_did_this_with %s"\n' % self.genome_version
        self.samfile.write(header)

    def add(self, chrom, strand, pos, cigar):
        """write a sam record, all other fields set to NA. pos is dna0.
        strand is 0 for pos and 1 for neg, or else just one of '+' or '-'.

        assumes file is open, use within 'with' block
        """

        # flags, ignore all the mate pair stuff and unmapped/duplicate indicators etc,
        # just use the strand field for now
        flags = 0x0000
        if strand == '-':
            flags = flags | 0x0010

        # construct record as tab delim string
        rec = ['na'] * 12
        rec[1] = str(flags)

        chrom = str(chrom)
        if chrom.startswith('chr'):
            chrom = chrom[3:]

        rec[2] = chrom
        rec[3] = str(pos + 1)  # sam wants dna1
        rec[5] = cigar

        self.samfile.write('\t'.join(rec) + '\n')


class SamBuilder(object):
    def __init__(self, genome_version, junctions, reads_dict=None):
        self.genome_version = genome_version
        self.junctions = junctions
        self.distribution = {}
        # this was str records but now tuples, redundant with reads_dict, should combine
        self.sam_records = []

        # keep track of all reads produced in given dict
        if reads_dict is not None:
            self.reads = reads_dict
        else:
            self.reads = {}

    def sim_reads(self, reads_per_junction, min_read_length, max_read_length, seed=None):

        if seed is not None:
            random.seed(seed)

        for junction in self.junctions:
            for _ in range(reads_per_junction):

                assert len(junction) > 0, "zero length junction! %s" % junction

                # pick a strand, 0=forward, 1=reverse
                strand = randint(0, 1)
                strand = '+-' [strand]

                # pick a read length
                l = randint(min_read_length, max_read_length)

                # pick a position within the read where the junction falls,
                # this forces at least one base on either side
                p = randint(1, l - 1)

                # compute left & right overhangs
                overhang_l, overhang_r = p, l - p

                # not sure if we should cap the soft-clip size,
                # for now we just let it be anywhere from 0 to 1-overhang size,

                sl = randint(0, overhang_l - 1)
                ml = overhang_l - sl

                sr = randint(0, overhang_r - 1)
                mr = overhang_r - sr

                assert ml > 0
                assert mr > 0
                assert l == sl + ml + mr + sr

                self.sam_records.append(self.sam_record(junction, strand, sl, ml, mr, sr))
                self.count_reads(junction, strand, sl, ml, mr, sr)

    def count_reads(self, junction, strand, sl, ml, mr, sr):
        """keep track of reads with a tuple of (chrom, strand, start, end, overhang_left, overhang_right)
        as dict keys ad count as value.  use this for later correctness accounting
        """

        key = (junction.chrom, strand, junction.start, junction.end, ml, mr)
        if key in self.reads:
            self.reads[key] += 1
        else:
            self.reads[key] = 1

    def write_sam(self, filename):
        with SamFile(filename, self.genome_version) as sf:
            for record in self.sam_records:
                sf.add(*record)

    def sam_record(self, junction, strand, sl, ml, mr, sr):
        """12 tab-delimited columns, we are only interested in second (flag), third (chrom),
        fourth (position), and sixth (cigar). we won't bother adding in sequence etc as that
        isn't used here.
        """

        # 1-junction cigar string, soft-clipping included only if non-zero
        if sl > 0:
            cigar = '%dS' % sl
        else:
            cigar = ''

        cigar = cigar + '%dM%dN%dM' % (ml, len(junction), mr)

        if sr > 0:
            cigar = cigar + '%dS' % sr

        pos = junction.expand(ml).start  # start of left-overhang in dna0
        return (junction.chrom, strand, pos, cigar)


sam_header1 = """
@PG	CL:"i_did_this_with hg19"
"""


def make_ralign(outfile, infiles, reference_genome=MiniGenome('hg19'), **kwargs):
    ReadAlignments.build_ralign(outfile, infiles, reference_genome, **kwargs)
    return reference_genome


class ReadAlignmentTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        fd, cls.tmpsam = mkstemp()
        os.close(fd)
        fd, cls.tmpralign = mkstemp()
        os.close(fd)

    @classmethod
    def tearDownClass(cls):  # Remove temporary files
        os.unlink(cls.tmpsam)
        os.unlink(cls.tmpralign)

    def make_ralign(self, reference_genome=Genome('hg19'), **kwargs):
        return make_ralign(self.tmpralign, [self.tmpsam], reference_genome, **kwargs)

    def test_invalid_init(self):
        hg19 = MiniGenome()
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            Junction('chr1', '+', 0, 0, hg19)
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            Alignment('chr1', '+', 0, 0, hg19)
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            AlignmentMatch('chr1', '+', 0, 0, hg19)

    def test_close(self):
        def test_access(self, child):
            with self.assertRaisesRegex(OSError, 'close'):
                child.chromosome

        ralign_file = "tests/data/read_alignments/single_variant.ralign"
        with ReadAlignments(ralign_file) as ralign:
            junctions = ralign.junctions
            alignments = ralign.alignments
            matches = ralign.matches
            variants = ralign.variants
            junction = junctions[0]
            alignment = alignments[0]
            match = matches[0]
            variant = variants[0]
            refg = variant.refg

        def test_find(self, table):
            interval = Interval("chr1", "+", 0, 0, refg)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_5p_aligned(interval)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_3p_aligned(interval)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_5p_within(interval)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_3p_within(interval)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_within(interval)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_overlapping(interval)
            with self.assertRaisesRegex(OSError, 'close'):
                table.find_exact(interval)

        test_find(self, junctions)
        test_find(self, alignments)
        test_find(self, matches)
        test_find(self, variants)
        test_access(self, junction)
        test_access(self, alignment)
        test_access(self, match)
        test_access(self, variant)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_only_junction(self):
        # junctions can only be determined by flanking aligned regions, but some SAMs have extracted only the junction
        # information from reads
        dumptext(
            self.tmpsam, sam_header1, """
        SRR1068905.2392677	83	1	155162027	60	2N	=	155161933	-170	CACTC	@ADD=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
        """)
        with self.assertRaisesRegex(ValueError, 'No alignment match'):
            self.make_ralign()


    def test_single_variant(self):
        dumptext(
            self.tmpsam, sam_header1, """
        SRR1068905.2392677	83	1	155162027	60	1M2N3I1M	=	155161933	-170	CACTC	@ADD=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
        """)

        ralign_file = "tests/data/read_alignments/single_variant.ralign"
        if "CI" not in os.environ:  # avoid downloading full 2bit on CI
            make_ralign(ralign_file, [self.tmpsam], Genome('hg19'))

        with ReadAlignments(ralign_file) as ralign:
            self.assertFalse(ralign.junctions.stranded)
            self.assertFalse(ralign.alignments.stranded)
            self.assertFalse(ralign.matches.stranded)
            self.assertFalse(ralign.variants.stranded)
            self.assertEqual(1, len(ralign.junctions))
            self.assertEqual(1, len(ralign.alignments))
            self.assertEqual(2, len(ralign.matches))
            self.assertEqual(2, len(ralign.variants))

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def testOne(self):
        hg19 = Genome('hg19')
        make_ralign(self.tmpralign, ["./tests/data/read_alignments/test2.sam"], hg19)

        target_junc_itv = Interval("chr4", "+", 674377, 674377, hg19)
        target_variant1 = Variant("chr4", 674372, "A", "T", hg19)
        target_variant2 = Variant("chr4", 674335, "A", "C", hg19)

        with ReadAlignments(self.tmpralign) as ralign:
            junctions = ralign.junctions
            alignments = ralign.alignments

            total_reads = alignments.find_overlapping(target_junc_itv)
            body_reads = [r for r in total_reads if len(r.matches) == 1]

            num_variants1_bodyreads = 0
            num_variants2_bodyreads = 0

            for read in body_reads:
                body_match = read.matches[0]
                variants = body_match.variants

                for v in variants:
                    if v == target_variant1:
                        num_variants1_bodyreads += 1

                    if v == target_variant2:
                        num_variants2_bodyreads += 1

            self.assertEqual(num_variants1_bodyreads, 44)
            self.assertEqual(num_variants2_bodyreads, 21)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_build(self):
        dumptext(
            self.tmpsam, sam_header1, """
        SRR1068905.2392674	83	1	155162021	60	47M1549N27M1I1M	=	155161960	-137	AAGTCTCCTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACAGCATGGTGGACGGCTGGGAGCTTTAGGAGG	EEDCFFFHGGHHHEHE@JJJJJJHCJJJJIIIEJJGGIIJJIIHJJGEIIIIIIIHGGJJJJIFHHGHFFFFFCCC	AS:i:-11	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:45C0C29	YS:i:0	YT:Z:CP	XS:A:+	NH:i:1
        SRR1068905.2392754	163	1	155162025	60	50M745N24M2D2M	=	155162849	155	CTCCTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACCTGTAACAACTCTTCCTCCTCCTCCCCAACCATTC	CCCFFFFFHHHHHIJJJJJJGGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIJJJJJJJJJJHHHFFFFEF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:42C33	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:42|S|rs4072037
        SRR1068905.2392676	83	1	155162027	60	41M1549N35M	=	155161837	-266	CCTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACAGCATGGTGGACGGCTGGGAGCTTTAGGAGGGGGCAC	?CECAC??B?;HEFJJJIGIGFGGGGIIGIGDEGIGEFGJJIGJIHGGJIJHF@IIJIHGHGGHDHFHDFDFD@@@	AS:i:-12	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:39C0C35	YS:i:0	YT:Z:CP	XS:A:+	NH:i:1
        SRR1068905.2392677	83	1	155162027	60	76M	=	155161933	-170	ACTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACCCGTAACAACTGTTGCGGGTTTAGGGGCTGTGGTAGC	@ADDDAC?D?7DDDDDDEIIIEDECDDIEDDBEDDDDDIEEIEDAECEIIIEEIIIIIEEEIIDDDDDDDA=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
        """)

        # 4th read a body read with 1 SN variant (155162027:C>A)

        hg19 = Genome('hg19')

        # -1 since sam is 1 based, and Interval is 0-based
        junctions_expected = [
            Interval("chr1", '+', 155162021 + 47 - 1, 155162021 + 47 + 1549 - 1, hg19),
            Interval("chr1", '+', 155162025 + 50 - 1, 155162025 + 50 + 745 - 1, hg19),
            # Interval("chr1", '+', 155162027 + 41 - 1, 155162027 + 41 + 1549 - 1, hg19) # same as element 0
        ]

        # alignment : start = col[3]-1  end = start + {all of M, D, N}
        alignments_expected = [
            Interval("chr1", '+', 155162021 - 1, 155162021 + 47 + 1549 + 27 + 0 + 1 - 1, hg19),
            Interval("chr1", '+', 155162025 - 1, 155162025 + 50 + 745 + 24 + 2 + 2 - 1, hg19),
            Interval("chr1", '+', 155162027 - 1, 155162027 + 41 + 1549 + 35 - 1, hg19),
            # A body read, cannot be accessed via junction
            Interval("chr1", '+', 155162027 - 1, 155162027 + 76 - 1, hg19),
        ]

        matches_expected = [
            Interval("chr1", '+', 155162021             - 1, 155162021 + 47                     - 1, hg19),
            Interval("chr1", '+', 155162021 + 47 + 1549 - 1, 155162021 + 47 + 1549 + 27 + 0 + 1 - 1, hg19),
            Interval("chr1", '+', 155162025             - 1, 155162025 + 50                     - 1, hg19),
            Interval("chr1", '+', 155162025 + 50 + 745  - 1, 155162025 + 50 + 745 + 24 + 2 + 2  - 1, hg19),
            Interval("chr1", '+', 155162027             - 1, 155162027 + 41                     - 1, hg19),
            Interval("chr1", '+', 155162027 + 41 + 1549 - 1, 155162027 + 41 + 1549 + 35         - 1, hg19),
            # A body read, entirety is a match
            Interval("chr1", '+', 155162027             - 1, 155162027 + 76                     - 1, hg19),
        ]  # yapf: disable

        variants_expected = [
            Variant.from_string('chr1:155162066:C:A', hg19),
            Variant.from_string('chr1:155162067:C:G', hg19),
            Variant.from_string('chr1:155163644::G', hg19),
            Variant.from_string('chr1:155162067:C:T', hg19),
            Variant.from_string('chr1:155162844:TC:', hg19),
            Variant.from_string('chr1:155162846:C:T', hg19),
            Variant.from_string('chr1:155162847:T:C', hg19),
            Variant.from_string('chr1:155162027:C:A', hg19),
        ]

        self.make_ralign(hg19)
        self.assertTrue(os.path.isfile(self.tmpralign), 'failed to create expected output ralign file')

        with ReadAlignments(self.tmpralign) as ralign:

            junctions = ralign.junctions
            alignments = ralign.alignments
            matches = ralign.matches
            variants = ralign.variants

            self.assertEqual(len(junctions), len(junctions_expected))  # =2, read 1, 3 share the same junction
            self.assertEqual(len(alignments), len(alignments_expected))  # =4, tmpsam has 4 rows, i.e. read alignments
            self.assertEqual(len(matches), len(matches_expected))  # =7, 4 reads, first 3 has 2 matches
            self.assertEqual(len(variants), len(variants_expected))  # =8

            # checking intervals are as expected
            self.assertEqual([junction.interval for junction in junctions], junctions_expected)
            self.assertEqual([alignment.interval for alignment in alignments], alignments_expected)
            self.assertEqual([match.interval for match in matches], matches_expected)
            self.assertEqual(list(variants), variants_expected)

            # check if each parent/child relationship is correctly mapped
            self.assertEqual(list(map(len, junctions)), list(map(len, junctions_expected)))
            self.assertEqual(list(map(len, alignments)), list(map(len, alignments_expected)))
            self.assertEqual(list(map(len, matches)), list(map(len, matches_expected)))

            self.assertEqual([len(junction.alignments) for junction in junctions], [2, 1])
            self.assertEqual([[alignment.interval for alignment in junction.alignments] for junction in junctions],
                             [[alignments_expected[0], alignments_expected[2]],
                              [alignments_expected[1]]])  # yapf: disable
            self.assertEqual([[match.interval for match in alignment.matches] for alignment in alignments],
                             [[matches_expected[0], matches_expected[1]],
                              [matches_expected[2], matches_expected[3]],
                              [matches_expected[4], matches_expected[5]],
                              [matches_expected[6]]])  # yapf: disable
            self.assertEqual([match.variants for match in matches],
                             [(variants_expected[0], variants_expected[1]),
                              (variants_expected[2], ),
                              (variants_expected[3], ),
                              (variants_expected[4], variants_expected[5], variants_expected[6]),
                              (variants_expected[0], variants_expected[1]),
                              (),
                              (variants_expected[7], )])  # yapf: disable

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_pileup(self):
        # identical to test_build + 1 N region on 5th row
        dumptext(
            self.tmpsam, sam_header1, """
        SRR1068905.2392674	83	1	155162021	60	47M1549N27M1I1M	=	155161960	-137	AAGTCTCCTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACAGCATGGTGGACGGCTGGGAGCTTTAGGAGG	EEDCFFFHGGHHHEHE@JJJJJJHCJJJJIIIEJJGGIIJJIIHJJGEIIIIIIIHGGJJJJIFHHGHFFFFFCCC	AS:i:-11	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:45C0C29	YS:i:0	YT:Z:CP	XS:A:+	NH:i:1
        SRR1068905.2392754	163	1	155162025	60	50M745N24M2D2M	=	155162849	155	CTCCTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACCTGTAACAACTCTTCCTCCTCCTCCCCAACCATTC	CCCFFFFFHHHHHIJJJJJJGGIIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIIJJJJJJJJJJHHHFFFFEF	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:42C33	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:42|S|rs4072037
        SRR1068905.2392676	83	1	155162027	60	41M1549N35M	=	155161837	-266	CCTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACAGCATGGTGGACGGCTGGGAGCTTTAGGAGGGGGCAC	?CECAC??B?;HEFJJJIGIGFGGGGIIGIGDEGIGEFGJJIGJIHGGJIJHF@IIJIHGHGGHDHFHDFDFD@@@	AS:i:-12	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:39C0C35	YS:i:0	YT:Z:CP	XS:A:+	NH:i:1
        SRR1068905.2392677	83	1	155162027	60	76M	=	155161933	-170	ACTTTTCTCCACCTGGGGTAGAGCTTGCATGACCAGAACCCGTAACAACTGTTGCGGGTTTAGGGGCTGTGGTAGC	@ADDDAC?D?7DDDDDDEIIIEDECDDIEDDBEDDDDDIEEIEDAECEIIIEEIIIIIEEEIIDDDDDDDA=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
        RefN	83	1	1	60	1M	=	0	0	A	@
        AltN	83	1	30000001	60	1M	=	0	0	N	@
        SRR1075603.69096319	99	6	31238879	60	75M	=	31238934	131	TCCAGGTGTCTGCGGAGCCACTCCACGCACAGGCCCTCCAGGTAGGCTCTCAGCTGCTCCGCCTCACGGGCCGCC	CCCFFFEFHHHHHJJJJJJJJJJJJJJJJJJIJJJJIJJJJJJJJJJJJJJHHHHHHFFFFDDDDDDDDDDDDDD	AS:i:-6	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:7A22G0T31G11	YS:i:-9	YT:Z:CP	NH:i:1	Zs:Z:7|S|rs77935220,23|S|rs1050685.1,31|S|rs2308590.1
        SRR1075603.69744482	163	6	31238237	60	26M587N50M	=	31238865	116	GTGGGTCACATGTGTCTTTGGGGGGTCCGCGCGCTGCAGCGTCTCCTTCCCGTTCTCCAGGTGTCTGCGGAGCCAC	@@	@FFFDDHHFDDEHIJJJJJFIJDDDDBDDDDDDDDDDCDBDBCDDCDDCDDBDDDDACDACDCCCDDDBDDBDD	AS:i:-10	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:9G14T2T34A13	YS:i:	0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:27|S|rs2308604,34|S|rs77935220
        """)

        hg19 = Genome('hg19')
        self.make_ralign(hg19)
        self.assertTrue(os.path.isfile(self.tmpralign), 'failed to create expected output ralign file')

        with ReadAlignments(self.tmpralign) as ralign:
            no_matches = Interval('chr1', '+', 155162027 + 76 - 1, 155162027 + 76 - 1, hg19).expand(0, 100)
            self.assertEqual([], ralign.matches.find_overlapping(no_matches))
            np.testing.assert_equal(np.zeros((len(no_matches), 5), np.int16), ralign.pileup(no_matches))

            deletion = Variant.from_string('chr1:155162844:TC:', hg19)
            self.assertEqual([deletion], ralign.variants.find_overlapping(deletion))
            np.testing.assert_equal(2 * [[0, 0, 0, 0, 1]], ralign.pileup(deletion))

            insertion = Variant.from_string('chr1:155163644::G', hg19)
            # need to look downstream since it's a repeat of Gs
            self.assertEqual([insertion], ralign.variants.find_overlapping(insertion.expand(1, 6)))
            np.testing.assert_equal([[0, 0, 1, 0, 0], [0, 1, 0, 0, 0]], ralign.pileup((insertion.expand(1, 6)))[5:])

            reference_nosample = Interval('chr1', '+', 155162027 + 41 + 1549 + 35 - 1, 155162027 + 41 + 1549 + 35 - 1,
                                          hg19).shift(-1).expand(0, 2)
            self.assertEqual(1, len(ralign.matches.find_overlapping(reference_nosample)))
            np.testing.assert_equal([[0, 1, 0, 0, 0], 5 * [0]], ralign.pileup(reference_nosample))

            variant_reference = Variant.from_string('chr1:155162027:C:A', hg19).expand(0, 1)
            self.assertEqual(4, len(ralign.matches.find_overlapping(variant_reference)))
            self.assertEqual(1, len(ralign.variants.find_overlapping(variant_reference)))
            np.testing.assert_equal([[1, 3, 0, 0, 0], [0, 4, 0, 0, 0]], ralign.pileup(variant_reference))

            overlapping_variants = Variant.from_string('chr1:155162067:C:G', hg19)
            self.assertEqual(4, len(ralign.matches.find_overlapping(overlapping_variants)))
            self.assertEqual(2, len(ralign.variants.find_overlapping(overlapping_variants)))
            np.testing.assert_equal([[0, 1, 2, 1, 0]], ralign.pileup(overlapping_variants))

            n_region = Interval('chr1', '+', 0, 2, hg19)
            self.assertEqual(1, len(ralign.matches.find_overlapping(n_region)))
            np.testing.assert_equal([[1, 0, 0, 0, 0], [0, 0, 0, 0, 0]], ralign.pileup(n_region))

            alt_n_region = Interval('chr1', '+', 30000000, 30000001, hg19)
            self.assertEqual(1, len(ralign.matches.find_overlapping(alt_n_region)))
            np.testing.assert_equal([[0, 0, 0, 0, 0]], ralign.pileup(alt_n_region))

            identical_variant = Variant.from_string('chr6:31238886:A:G', hg19)
            self.assertEqual(2, len(ralign.alignments.find_overlapping(identical_variant)))
            self.assertEqual(2, len(ralign.matches.find_overlapping(identical_variant)))
            self.assertEqual(1, len(ralign.variants.find_overlapping(identical_variant)))
            np.testing.assert_equal([[0, 0, 2, 0, 0]], ralign.pileup(identical_variant))

    def test_pileup_indices(self):
        np.testing.assert_equal([0, 1, 2, 3, 0, 1, 2, 3, 255], ReadAlignments._dna_as_pileup_index('acgtACGTN'))

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_duplicates(self):
        dumptext(
            self.tmpsam, sam_header1, """
        SRR1068905.2392677	1024	1	155162027	60	1M2N3I1M	=	155161933	-170	CACTC	@ADD=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
        """)

        self.make_ralign()
        with ReadAlignments(self.tmpralign) as ralign:
            self.assertEqual(0, len(ralign.alignments))

        self.make_ralign(include_duplicates=True)
        with ReadAlignments(self.tmpralign) as ralign:
            self.assertEqual(1, len(ralign.alignments))


def make_jralign(outfile, infiles, reference_genome=MiniGenome('hg19'), **kwargs):
    JReadAlignments.build_jralign(outfile, infiles, reference_genome, **kwargs)
    return reference_genome


class JReadAlignmentTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

        fd, cls.tmpsam = mkstemp()
        os.close(fd)
        fd, cls.tmpralign = mkstemp()
        os.close(fd)

    @classmethod
    def tearDownClass(cls):  # Remove temporary files
        os.unlink(cls.tmpsam)
        os.unlink(cls.tmpralign)

    def make_jralign(self, reference_genome=MiniGenome('hg19'), **kwargs):
        return make_jralign(self.tmpralign, [self.tmpsam], reference_genome, **kwargs)

    def test_only_junction(self):
        # junctions can only be determined by flanking aligned regions, but some SAMs have extracted only the junction
        # information from reads
        dumptext(
            self.tmpsam, sam_header1, """
            SRR1068905.2392677	83	1	155162027	60	2N	=	155161933	-170	CACTC	@ADD=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
            """)

        self.make_jralign()
        with JReadAlignments(self.tmpralign) as jralign:
            self.assertFalse(jralign.junctions.stranded)
            self.assertEqual(1, len(jralign.junctions))
            self.assertEqual(1, len(jralign.junctions[0]))
            self.assertEqual(0, jralign.junctions[0][0].left)
            self.assertEqual(0, jralign.junctions[0][0].right)
            self.assertEqual(0, jralign.junctions[0][0].num_variants)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_variant_feature(self):
        genome = Genome("hg19")

        def find_variants(chrom, pos, cigar, seq, genome):
            """ Finds variants based on CIGAR string and SEQ field in SAM """

            consume_seq = ['M', 'I', 'S']
            consume_ref = ['M', 'D', 'N']

            ref_pos, seq_pos = pos, 0
            variants = set()

            for code_len, code in re.findall("([0-9]+)(M|I|S|D|N){1}", cigar):
                code_len = int(code_len)
                if code == 'M':
                    itv = Interval(chrom, "+", ref_pos, ref_pos + code_len, genome)
                    dna = genome.dna(itv.shift(-1))

                    for i, v in enumerate(dna):
                        ref_nt = dna[i]
                        seq_nt = seq[seq_pos + i]
                        if ref_nt != seq_nt:
                            # print("ref_dna = {} [{}-{})".format(dna, ref_pos, ref_pos+code_len))
                            # print("seq_dna = {} [{}-{})".format(seq[seq_pos:seq_pos+code_len], seq_pos, seq_pos+code_len))
                            # print("found variant {} > {} ({})".format(ref_nt, seq_nt, ref_pos+i), file=sys.stderr)
                            v = Variant(chrom, ref_pos + i - 1, ref_nt, seq_nt, genome)
                            variants.add(v)

                if code in consume_seq:
                    seq_pos += code_len
                if code in consume_ref:
                    ref_pos += code_len

            return variants

        sam_file = os.path.join(TEST_DATA_FOLDER, 'test.sam')

        make_jralign(self.tmpralign, [sam_file], genome, include_variants=True)
        self.assertTrue(os.path.isfile(self.tmpralign), 'failed to create expected output jralign file')

        # Collects variants for each read alignments
        expected_variants = set()
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'): continue
                split_line = line.split('\t')
                chrom, pos, cigar, seq = split_line[2], int(split_line[3]), split_line[5], split_line[9]
                expected_variants.update(find_variants(chrom, pos, cigar, seq, genome))

        # tests
        with JReadAlignments(self.tmpralign) as jralign:

            junc_table = jralign.junctions
            self.assertIsInstance(junc_table, JunctionReadAlignmentsTable)
            variant_table = jralign.variants
            self.assertIsInstance(variant_table, VariantTable)

            # Check variant table matches expected set of variants in its entirety
            self.assertEqual(set(variant_table), expected_variants)

            for junc in junc_table:
                self.assertIsInstance(junc, JunctionReadAlignments)

                for read in junc:
                    self.assertIsInstance(read, JunctionReadAlignment)
                    read_variants = read.variants()
                    self.assertEqual(read.num_variants, len(read_variants))

                    # Check that variants pulled from read alignments are also in the expected variants
                    for v in read_variants:
                        self.assertIn(v, expected_variants)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_contingency_matrix_calculation(self):
        hg19 = Genome('hg19')

        sam_file = os.path.join(TEST_DATA_FOLDER, 'test.sam')

        make_jralign(self.tmpralign, [sam_file], hg19, include_variants=True)

        reads_const_ref = 0
        reads_const_alt = 0
        reads_alter_ref = 0
        reads_alter_alt = 0

        the_var = Variant.from_string("chr1:155162067:C:T", hg19)
        const_interval = Interval("chr1", '+', 155162074, 155162076, hg19)
        alter_interval = Interval("chr1", '+', 155162101, 155162103, hg19)

        with JReadAlignments(self.tmpralign) as jralign:
            const_juncs = jralign.junctions.find_5p_within(const_interval)
            alter_juncs = jralign.junctions.find_5p_within(alter_interval)

            # for reads supporting the constitutive junction
            reads_on_var = 0
            for const_junc in const_juncs:
                for read in const_junc:
                    if the_var.position > const_junc.start - read.left:
                        reads_on_var += 1
                    for var in read.variants():
                        if (var.start == the_var.position and var.ref == the_var.ref and var.alt == the_var.alt):
                            reads_const_alt += 1
            reads_const_ref = reads_on_var - reads_const_alt

            # for reads supporting alternative junction
            reads_on_var = 0
            for alter_junc in alter_juncs:
                for read in alter_junc:
                    if the_var.position > alter_junc.start - read.left:
                        reads_on_var += 1
                    for var in read.variants():
                        if (var.start == the_var.position and var.ref == the_var.ref and var.alt == the_var.alt):
                            reads_alter_alt += 1
            reads_alter_ref = reads_on_var - reads_alter_alt

            contingency_2x2_table = np.array([[reads_const_ref, reads_const_alt],
                                                [reads_alter_ref, reads_alter_alt]])
            self.assertTrue(np.array_equal(np.array([[0, 26], [8, 1]]), contingency_2x2_table))

    def test_close(self):
        num_junctions = 1
        num_reads_per_junction = 1
        min_read_length = 50
        max_read_length = 100

        junctions = sample_junctions(num_junctions)

        builder = SamBuilder('hg19', junctions)
        builder.sim_reads(num_reads_per_junction, min_read_length, max_read_length)

        # write out the sam file
        builder.write_sam(self.tmpsam)

        # build jralign
        genome = self.make_jralign()
        with JReadAlignments(self.tmpralign) as jralign:
            junctions = jralign.junctions
            junction = junctions[0]

        interval = Interval("chr1", "+", 0, 0, genome)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_5p_aligned(interval)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_3p_aligned(interval)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_5p_within(interval)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_3p_within(interval)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_within(interval)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_overlapping(interval)
        with self.assertRaisesRegex(OSError, 'close'):
            junctions.find_exact(interval)

        with self.assertRaisesRegex(OSError, 'close'):
            junction.chromosome

    def test_one_sam(self):
        """detailed correctness test, builds an jralign from a programatically
        generated sam file, make sure all reads are accounted for and correct
        """

        num_junctions = 10
        num_reads_per_junction = 100
        min_read_length = 50
        max_read_length = 100

        junctions = sample_junctions(num_junctions)

        builder = SamBuilder('hg19', junctions)
        builder.sim_reads(num_reads_per_junction, min_read_length, max_read_length)

        # write out the sam file
        builder.write_sam(self.tmpsam)

        # build jralign
        self.make_jralign()

        self.assertTrue(os.path.isfile(self.tmpralign), 'failed to create expected output jralign file')

        # load & validate the jralign
        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(num_junctions, len(jralign.junctions))

            reads = builder.reads

            # iterate over all the reads in and check
            for alignment in jralign.junctions:
                self.assertEqual(num_reads_per_junction, alignment.num_reads)
                assert alignment.strand == '+'  # this might change in future

                for read in alignment:
                    # this key gets matched by the one produced in the builder, validates
                    # all positions related data
                    key = (alignment.chrom, read.strand, alignment.start, alignment.end, read.left, read.right)

                    # read must be in our accounting dict, verify & remove
                    self.assertIn(key, reads, 'alignment file contains unexpected read: %s' % str(key))
                    count = reads[key]
                    if count == 1:
                        del reads[key]
                    else:
                        reads[key] = count - 1

            # we should have completely drained the accounting dict, verify
            self.assertEqual(0, len(reads), 'alignment file is missing reads')

    def test_multi_sam(self):
        """handful of small sam files. each different. but likely no overlapping junctions
        """
        num_sam_files = 5
        num_junctions = 10
        num_reads_per_junction = 10
        min_read_length = 50
        max_read_length = 100

        reads = {}  # aggregates all simulated reads

        with TempDir() as tmpdir:
            sam_files = []
            for sammy in range(num_sam_files):
                junctions = sample_junctions(num_junctions)
                builder = SamBuilder('hg19', junctions, reads)
                builder.sim_reads(num_reads_per_junction, min_read_length, max_read_length)
                sam_file = tmpdir('simulated-%d.sam' % sammy)
                sam_files.append(sam_file)
                builder.write_sam(sam_file)

            # build jralign
            make_jralign(self.tmpralign, sam_files)

            self.assertTrue(os.path.isfile(self.tmpralign), 'failed to create expected output jralign file')

            # load & validate the jralign
            '''
            with JReadAlignments(jralign_file) as jralign:
                # iterate over all the reads in and check
                for alignment in jralign.junctions:
                    self.assertEqual(num_reads_per_junction, alignment.num_reads)
                    assert alignment.strand == '+'  # this might change in future

                    for read in alignment:
                        # this key gets matched by the one produced in the builder, validates
                        # all positions related data
                        key = (alignment.chrom, read.strand, alignment.start, alignment.end,
                               read.left, read.right)

                        # read must be in our accounting dict, verify & remove
                        self.assertIn(key, reads, 'jralign contains unexpected read: %s' % str(key))
                        count = reads[key]
                        if count == 1:
                            del reads[key]
                        else:
                            reads[key] = count-1

                # we should have completely drained the accounting dict, verify
                jralign.close()
                self.assertEqual(0, len(reads), 'alignment file is missing reads')
            '''

    def test_duplicates(self):
        dumptext(
            self.tmpsam, sam_header1, """
        SRR1068905.2392677	1024	1	155162027	60	1M2N3I1M	=	155161933	-170	CACTC	@ADD=A:1+	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:40C35	YS:i:0	YT:Z:CP	XS:A:-	NH:i:1	Zs:Z:40|S|rs4072037
        """)

        self.make_jralign()
        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(0, len(jralign.junctions))

        self.make_jralign(include_duplicates=True)
        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))


    def test_allow(self):
        sam_file = os.path.join(TEST_DATA_FOLDER, 'exclude.sam')

        # no allow
        genome = make_jralign(self.tmpralign, [sam_file], exclude=[], allow=[])

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(3, len(jralign.junctions))

        # allow stretching completely over one junction,
        allow = [Interval('chr19', '+', 35551529 - 1000, 35551529 + 1000, genome)]
        # allow = [Interval('chr19', '+', 35551417, 35551417, 'genome)]
        make_jralign(self.tmpralign, [sam_file], genome, allow=allow)

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))

        # same allow as above, with also a exclude of a
        # 0-len interval within junction, junction is blocked
        allow = [Interval('chr19', '+', 35551529 - 1000, 35551529 + 1000, genome)]
        exclude = [Interval('chr19', '+', 35551417, 35551417, genome)]
        make_jralign(self.tmpralign, [sam_file], genome, exclude=exclude, allow=allow)

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(0, len(jralign.junctions))

    # Note: for all the insertion/deletion tests,
    # our rules is that overhang calculations include
    # the insertions but not deletions, whereas
    # junction positions include the deletions but not
    # insertions. That is, junctions positions are in
    # reference coordinates, and overhang sizes are
    # private to the sample genome
    #
    # so if we have 1M1D1M2N1M (allowing 1 base matches
    # for illustration), at position 2, we have:
    #
    #    0 1 2 3 4 5 6 7 8 9
    #
    #        M D M - - M
    #
    # so the start coordinate of the junction is 2+1+1+1 = 5,
    # that is the D at position 3 is counted, but the left overhang
    # is only 2 with the D not being counted.
    #
    # insertions work similarly, so with 1M1I1M2N1M at pos 2:
    #
    #    0 1 2 3 4 5 6 7 8 9
    #
    #        M M - - M
    #         ^
    #         I
    #
    #  here the junction has coordinates 2+1+1 = 4 (the insertion squeezed
    # in the diagram between the M's isn't counted), but the left
    # overhang is 3 with the insertion being counted
    #
    def test_insertion_left(self):
        with SamFile(self.tmpsam, 'hg19') as sf:
            sf.add('chr1', '+', 1000, '1I1M10N2M')

        self.make_jralign()

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))

            junction = jralign.junctions[0]
            self.assertEqual('chr1', junction.chrom)
            self.assertEqual(1000 + 1, junction.start)
            self.assertEqual(1000 + 1 + 10, junction.end)

            self.assertEqual(1, len(junction))

            read = junction[0]
            self.assertEqual('+', read.strand)
            self.assertEqual(1, read.left)
            self.assertEqual(2, read.right)

    def test_insertion_right(self):
        with SamFile(self.tmpsam, 'hg19') as sf:
            sf.add('chr1', '+', 1000, '1M10N1M1I')

        self.make_jralign()

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))

            junction = jralign.junctions[0]
            self.assertEqual('chr1', junction.chrom)
            self.assertEqual(1000 + 1, junction.start)
            self.assertEqual(1000 + 1 + 10, junction.end)

            self.assertEqual(1, len(junction))

            read = junction[0]
            self.assertEqual('+', read.strand)
            self.assertEqual(1, read.left)
            self.assertEqual(1, read.right)

    def test_deletion_left(self):
        with SamFile(self.tmpsam, 'hg19') as sf:
            sf.add('chr1', '+', 1000, '1D1M10N1M')

        self.make_jralign()

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))

            junction = jralign.junctions[0]
            self.assertEqual('chr1', junction.chrom)
            self.assertEqual(1000 + 1 + 1, junction.start)
            self.assertEqual(1000 + 1 + 1 + 10, junction.end)

            self.assertEqual(1, len(junction))

            read = junction[0]
            self.assertEqual('+', read.strand)
            self.assertEqual(1, read.left)  # only 1, the 1D doesn't count
            self.assertEqual(1, read.right)

    def test_deletion_right(self):
        with SamFile(self.tmpsam, 'hg19') as sf:
            sf.add('chr1', '+', 1000, '1M10N1M1D')

        self.make_jralign()

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))

            junction = jralign.junctions[0]
            self.assertEqual('chr1', junction.chrom)
            self.assertEqual(1000 + 1, junction.start)
            self.assertEqual(1000 + 1 + 10, junction.end)

            self.assertEqual(1, len(junction))

            read = junction[0]
            self.assertEqual('+', read.strand)
            self.assertEqual(1, read.left)
            self.assertEqual(1, read.right)  # only 1, the 1D doesn't count

    def test_multibase_deletion(self):
        """deletion of 2 bases
        """
        with SamFile(self.tmpsam, 'hg19') as sf:
            sf.add('chr1', '+', 1000, '2D1M10N1M')

        self.make_jralign()

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))

            junction = jralign.junctions[0]
            self.assertEqual('chr1', junction.chrom)
            self.assertEqual(1000 + 2 + 1, junction.start)
            self.assertEqual(1000 + 2 + 1 + 10, junction.end)

            self.assertEqual(1, len(junction))

            read = junction[0]
            self.assertEqual('+', read.strand)
            self.assertEqual(1, read.left)  # only 1, the 2D doesn't count
            self.assertEqual(1, read.right)

    def test_multiple_junctions(self):
        """two N sections in cigar string.

        reads spanning two junctions are split into a pair of reads,
        with the overhang of the central match being added to the
        total overhang in both cases (on opposite ends).
        """
        with SamFile(self.tmpsam, 'hg19') as sf:
            sf.add('chr1', '-', 1000, '1M10N2M10N3M')

        self.make_jralign()

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(2, len(jralign.junctions))  # expect two junctions

            # need to order these, make sure junction_l is the left-most
            junction_l = jralign.junctions[0]
            junction_r = jralign.junctions[1]

            if junction_l.start > junction_r.start:
                junction_l, junction_r = junction_r, junction_l

            # check first (left) junction
            self.assertEqual('chr1', junction_l.chrom)
            self.assertEqual(1000 + 1, junction_l.start)
            self.assertEqual(1000 + 1 + 10, junction_l.end)

            self.assertEqual(1, len(junction_l))  # 1 read

            read = junction_l[0]
            self.assertEqual('-', read.strand)
            self.assertEqual(1, read.left)
            self.assertEqual(2 + 3, read.right)  # middle + right

            # check second (right)
            self.assertEqual('chr1', junction_r.chrom)
            self.assertEqual(1000 + 1 + 10 + 2, junction_r.start)
            self.assertEqual(1000 + 1 + 10 + 2 + 10, junction_r.end)

            self.assertEqual(1, len(junction_r))  # 1 read

            read = junction_r[0]
            self.assertEqual('-', read.strand)
            self.assertEqual(1 + 2, read.left)  # left + middle
            self.assertEqual(3, read.right)

    def test_large_overhang(self):
        sam_file = os.path.join(TEST_DATA_FOLDER, 'large_overhang.sam')

        # overhang of 255 is ok
        make_jralign(self.tmpralign, [sam_file])

        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))
            self.assertEqual(1, len(jralign.junctions[0]))

            self.assertEqual(255, jralign.junctions[0][0].right)

        sam_file = os.path.join(TEST_DATA_FOLDER, 'too_large_overhang.sam')

        # overhang of 256 is an error, but handled gracefully with exception
        with self.assertRaises(ValueError):
            make_jralign(self.tmpralign, [sam_file])

        # clamping overhang to 255
        make_jralign(self.tmpralign, [sam_file], overhang_error="clamp")
        with JReadAlignments(self.tmpralign) as jralign:
            self.assertEqual(1, len(jralign.junctions))
            self.assertEqual(1, len(jralign.junctions[0]))

            self.assertEqual(255, jralign.junctions[0][0].right)


class SimpleReadDistTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

    def makedist(self, sam_file, jralign_file, rdist_file):
        """Helper, builds rdist

        sam_file pre-exists and is given as a full  path, jralign_file and
        rdist_file are the names of the new files to be created
        """

        # build jralign
        make_jralign(jralign_file, [sam_file])

        self.assertTrue(os.path.isfile(jralign_file), 'failed to create expected output jralign file')

        # build an rdist
        ReadDistributions.build_rdist(rdist_file, [jralign_file])

        self.assertTrue(os.path.isfile(rdist_file), 'failed to create expected output rdist file')

    def test_invalid_init(self):
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            JunctionReadAlignments('chr1', '+', 0, 0, MiniGenome())
        with self.assertRaisesRegex(RuntimeError, "internal type"):
            JunctionReadDistribution('chr1', '+', 0, 0, MiniGenome())

    def test_close(self):
        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, '1bp_overhang.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            self.makedist(sam_file, jralign_file, rdist_file)

            with ReadDistributions(rdist_file) as rdist:
                junctions = rdist.junctions
                junction = rdist.junctions[0]

            interval = Interval("chr1", "+", 0, 0, "hg19")
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_5p_aligned(interval)
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_3p_aligned(interval)
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_5p_within(interval)
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_3p_within(interval)
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_within(interval)
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_overlapping(interval)
            with self.assertRaisesRegex(IOError, 'close'):
                junctions.find_exact(interval)

            with self.assertRaisesRegex(IOError, 'close'):
                junction.num_reads()

    def test_1bp_overhang(self):

        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, '1bp_overhang.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            self.makedist(sam_file, jralign_file, rdist_file)

            with ReadDistributions(rdist_file) as rdist:
                self.assertFalse(rdist.junctions.stranded)

                # expect 1 junction, 2 readcount objects, each with count 1.
                # num_reads for the junction is two since we count reads at
                # each end of the junction separately
                self.assertEqual(1, len(rdist.junctions))

                junction = rdist.junctions[0]  # this gives # of count objects, same as num_counts
                self.assertEqual(2, len(junction))
                self.assertEqual(2, junction.num_counts)
                self.assertEqual(1, junction.num_reads)

                # jrcounts don't seem to be sortable (though they don't give an error either)
                # copy into list of (strand, shift, count) tuples, sort, and use for validation
                counts = [(c.strand, c.shift, c.count) for c in junction]
                counts.sort()

                self.assertEqual([('+', -1, 1), ('+', 1, 1)], counts)

    def test_2reads_opp_strands(self):

        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, '2reads_opp_strands.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            self.makedist(sam_file, jralign_file, rdist_file)

            with ReadDistributions(rdist_file) as rdist:

                # expect 1 junction, 2 readcount objects, each with count 1.
                # num_reads for the junction is two since we count reads at
                # each end of the junction separately
                self.assertEqual(1, len(rdist.junctions))

                junction = rdist.junctions[0]
                self.assertEqual(4, len(junction))  # this gives # of count objects, same as num_counts
                self.assertEqual(4, junction.num_counts)
                self.assertEqual(2, junction.num_reads)

                # jrcounts don't seem to be sortable (though they don't give an error either)
                # copy into list of (strand, shift, count) tuples, sort, and use for validation
                counts = [(c.strand, c.shift, c.count) for c in junction]
                counts.sort()

                self.assertEqual([('+', -1, 1), ('+', 2, 1), ('-', -3, 1), ('-', 4, 1)], counts)

    def test_shift_match(self):
        """couple of reads with same overhang, count > 1
        """
        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, '2reads_same_overhang.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            self.makedist(sam_file, jralign_file, rdist_file)

            with ReadDistributions(rdist_file) as rdist:

                # expect 1 junction, 2 readcount objects, each with count 1.
                # num_reads for the junction is two since we count reads at
                # each end of the junction separately
                self.assertEqual(1, len(rdist.junctions))

                junction = rdist.junctions[0]
                self.assertEqual(2, len(junction))  # this gives # of count objects, same as num_counts
                self.assertEqual(2, junction.num_counts)
                self.assertEqual(2, junction.num_reads)

                # jrcounts don't seem to be sortable (though they don't give an error either)
                # copy into list of (strand, shift, count) tuples, sort, and use for validation
                counts = [(c.strand, c.shift, c.count) for c in junction]
                counts.sort()

                # notice we are implicitly testing that the sum of the counts is
                # 4 (last elem of each tuple)
                self.assertEqual([('+', -1, 2), ('+', 2, 2)], counts)

    def test_min_overhang(self):
        """apply min_overhang threshold
        """

        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, 'min_overhang.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            make_jralign(jralign_file, [sam_file], min_overhang=2)

            ReadDistributions.build_rdist(rdist_file, [jralign_file])

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(1, len(rdist.junctions))

                junction = rdist.junctions[0]
                self.assertEqual(6, len(junction))
                self.assertEqual(6, junction.num_counts)
                self.assertEqual(3, junction.num_reads)

                counts = [(c.strand, c.shift, c.count) for c in junction]
                counts.sort()

                self.assertEqual([('+', -2, 1), ('+', 3, 1), ('-', -4, 1), ('-', -3, 1), ('-', 4, 1), ('-', 5, 1)],
                                 counts)
            # rebuild rdist from same jralign, this time with a threshold
            ReadDistributions.build_rdist(rdist_file, [jralign_file], min_overhang=3)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(1, len(rdist.junctions))

                junction = rdist.junctions[0]
                self.assertEqual(4, len(junction))
                self.assertEqual(4, junction.num_counts)
                self.assertEqual(2, junction.num_reads)

                counts = [(c.strand, c.shift, c.count) for c in junction]
                counts.sort()

                self.assertEqual([('-', -4, 1), ('-', -3, 1), ('-', 4, 1), ('-', 5, 1)], counts)

    def test_min_reads(self):
        """apply min_overhang threshold
        """

        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, 'min_reads.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            # jralign with a threshold
            make_jralign(jralign_file, [sam_file], min_reads=2)

            # rdist with no threshold
            ReadDistributions.build_rdist(rdist_file, [jralign_file])

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(2, len(rdist.junctions))

            # rebuild rdist from same jralign, this time with a threshold
            ReadDistributions.build_rdist(rdist_file, [jralign_file], min_reads=3)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(1, len(rdist.junctions))

    def test_pass_filter(self):
        """emulate min_reads with pass_filter
        """

        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, 'min_reads.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            make_jralign(jralign_file, [sam_file])
            with JReadAlignments(jralign_file) as jralign:
                def min_three_reads(infile, junction, read):
                    return jralign.junctions[junction].num_reads >= 3

                with ReadDistributions.from_jraligns([jralign_file],
                                                     min_reads=1,
                                                     pass_filter=min_three_reads,
                                                     outfile=rdist_file) as rdist:
                    self.assertEqual(1, len(rdist.junctions))

    def test_exclude(self):
        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, 'exclude.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            # no exclude
            genome = make_jralign(jralign_file, [sam_file], exclude=[])

            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(3, len(jralign.junctions))

            # exclude stretching completely over one junction
            exclude = [Interval('chr19', '+', 35551529 - 1000, 35551529 + 1000, genome)]
            make_jralign(jralign_file, [sam_file], exclude=exclude)
            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(2, len(jralign.junctions))

            # exclude interval of 0 size on the left edge of the junction,
            # does not block the junction
            exclude = [Interval('chr19', '+', 35551416, 35551416, genome)]
            make_jralign(jralign_file, [sam_file], exclude=exclude)

            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(3, len(jralign.junctions))

            # exclude interval of 0 size at first position within junction
            # does block interval
            exclude = [Interval('chr19', '+', 35551417, 35551417, genome)]
            make_jralign(jralign_file, [sam_file], exclude=exclude)

            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(2, len(jralign.junctions))

            # exclude interval of 0 size on right edge of junction
            # does not block interval
            exclude = [Interval('chr19', '+', 35551530, 35551530, genome)]
            make_jralign(jralign_file, [sam_file],exclude=exclude)

            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(3, len(jralign.junctions))

            # exclude interval of 0 size at rightmost position within junction,
            # does block interval
            exclude = [Interval('chr19', '+', 35551529, 35551529, genome)]
            make_jralign(jralign_file, [sam_file], exclude=exclude)

            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(2, len(jralign.junctions))

            # jralign with no exclude, test excludeing within rdist
            # first no exclude
            make_jralign(jralign_file, [sam_file], exclude=[])

            ReadDistributions.build_rdist(rdist_file, [jralign_file], exclude=[])

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(3, len(rdist.junctions))

            # exclude stretching completely over one junction
            exclude = [Interval('chr19', '+', 35551529 - 1000, 35551529 + 1000, genome)]
            ReadDistributions.build_rdist(rdist_file, [jralign_file], exclude=exclude)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(2, len(rdist.junctions))

            # exclude interval of 0 size on the left edge of the junction,
            # does not block the junction
            exclude = [Interval('chr19', '+', 35551416, 35551416, genome)]
            ReadDistributions.build_rdist(rdist_file, [jralign_file], exclude=exclude)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(3, len(rdist.junctions))

            # exclude interval of 0 size at first position within junction
            # does block interval
            exclude = [Interval('chr19', '+', 35551417, 35551417, genome)]
            ReadDistributions.build_rdist(rdist_file, [jralign_file], exclude=exclude)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(2, len(rdist.junctions))

            # exclude interval of 0 size on right edge of junction
            # does not block interval
            exclude = [Interval('chr19', '+', 35551530, 35551530, genome)]
            ReadDistributions.build_rdist(rdist_file, [jralign_file], exclude=exclude)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(3, len(rdist.junctions))

            # exclude interval of 0 size at rightmost position within junction,
            # does block interval
            exclude = [Interval('chr19', '+', 35551529, 35551529, genome)]

            ReadDistributions.build_rdist(rdist_file, [jralign_file], exclude=exclude)

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(2, len(rdist.junctions))

    def test_file_handle_regression(self):

        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, '1bp_overhang.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')

            self.makedist(sam_file, jralign_file, rdist_file)

            for trial in range(20000):
                with ReadDistributions(rdist_file):
                    pass  # Make sure we can open/close the file 20000 times with no too-many-files errors

    def test_table_indexing(self):
        with TempDir() as tmpdir:
            sam_file = os.path.join(TEST_DATA_FOLDER, '3junctions.sam')
            jralign_file = tmpdir('test.jralign')
            rdist_file = tmpdir('test.rdist')
            self.makedist(sam_file, jralign_file, rdist_file)

            # while we have an jralign table open, check that indexing of a set of read alignments is Pythonic
            with JReadAlignments(jralign_file) as jralign:
                check_pythonic_indexing(self, jralign.junctions)  # Test JunctionReadAlignmentsTable indexing
                check_pythonic_indexing(self, jralign.junctions[0])  # Test JunctionReadAlignments indexing

            # while we have an rdist table open, check that indexing of a read distribution is Pythonic
            with ReadDistributions(rdist_file) as rdist:
                check_pythonic_indexing(self, rdist.junctions)  # Test JunctionReadDistributionTable indexing
                check_pythonic_indexing(self, rdist.junctions[0])  # Test JunctionReadDistribution

    def test_big_reads(self):
        """test > 2**16 reads mapping to a junction. because rdist has counters for each side,
        maybe we go > 2**17
        """

        # pick a specific junction, determine its length, then programmatically
        # generate lots of reads with 1 overhang on each side, so they will
        # all accumulate in the dist
        junction = MiniGenome().introns[0]
        junctions = [junction]

        num_reads_per_junction = 2**17 + 10
        min_read_length = max_read_length = 2

        builder = SamBuilder('hg19', junctions)
        builder.sim_reads(num_reads_per_junction, min_read_length, max_read_length)

        with TempDir() as tmpdir:
            # write out the sam file
            sam_file = tmpdir('simulated.sam')
            builder.write_sam(sam_file)

            # build jralign
            jralign_file = tmpdir('test.jralign')
            make_jralign(jralign_file, [sam_file])

            self.assertTrue(os.path.isfile(jralign_file), 'failed to create expected output jralign file')

            # load & validate the jralign
            with JReadAlignments(jralign_file) as jralign:
                self.assertEqual(1, len(jralign.junctions))

                reads = builder.reads

                # iterate over all the reads in and check
                for alignment in jralign.junctions:
                    self.assertEqual(num_reads_per_junction, alignment.num_reads)
                    assert alignment.strand == '+'  # this might change in future

                    for read in alignment:
                        # this key gets matched by the one produced in the builder, validates
                        # all positions related data
                        key = (alignment.chrom, read.strand, alignment.start, alignment.end, read.left, read.right)

                        # read must be in our accounting dict, verify & remove
                        self.assertIn(key, reads, 'alignment file contains unexpected read: %s' % str(key))
                        count = reads[key]
                        if count == 1:
                            del reads[key]
                        else:
                            reads[key] = count - 1

                # we should have completely drained the accounting dict, verify
                self.assertEqual(0, len(reads), 'alignment file is missing reads')

            # build the rdist and verify
            rdist_file = tmpdir('test.rdist')
            ReadDistributions.build_rdist(rdist_file, [jralign_file])

            self.assertTrue(os.path.isfile(rdist_file), 'failed to create expected output jralign file')

            with ReadDistributions(rdist_file) as rdist:
                self.assertEqual(1, len(rdist.junctions))

                junction = rdist.junctions[0]
                self.assertEqual(4, len(junction))
                self.assertEqual(4, junction.num_counts)
                self.assertEqual(131082, junction.num_reads)  # == (2**17+10)

                # jrcounts don't seem to be sortable (though they don't give an error either)
                # copy into list of (strand, shift, count) tuples, sort, and use for validation
                counts = [(c.strand, c.shift) for c in junction]
                counts.sort()

                # counts here omits the actual count, check that separately since we can't assert
                # the individual values, just the sum
                self.assertEqual([('+', -1), ('+', 1), ('-', -1), ('-', 1)], counts)
                self.assertEqual(262164, sum(c.count for c in junction))


class StrandedTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "assertRaisesRegex"):
            cls.assertRaisesRegex = cls.assertRaisesRegexp

    def setUp(self):
        fd, self.tmpsam = tempfile.mkstemp(suffix='.sam')
        os.close(fd)
        fd, self.tmpout = tempfile.mkstemp()
        os.close(fd)

    def tearDown(self):
        os.unlink(self.tmpsam)
        os.unlink(self.tmpout)

    def make_ralign(self, reference_genome=MiniGenome('hg19'), **kwargs):
        return make_ralign(self.tmpout, [self.tmpsam], reference_genome, **kwargs)

    def make_jralign(self, reference_genome=MiniGenome('hg19'), **kwargs):
        return make_jralign(self.tmpout, [self.tmpsam], reference_genome, **kwargs)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full anno on CI")
    def test_detect_with_anno_miss(self):
        dumptext(
            self.tmpsam, sam_header1, """
        Miss	0	chr1	0	255	1M1N1M	*	0	0	*	*
        """)

        with self.assertRaisesRegex(ValueError, 'cannot align'):
            self.make_ralign(Genome('gencode.v19'))

    @unittest.skipIf('CI' in os.environ, "avoid downloading full anno on CI")
    def test_detect_with_anno_ambiguous(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ENST00000466557.2.Intron6	0	chr1	155831	255	1M8431N1M	*	0	0	*	*
        ENST00000496488.1.Intron1	0	chr1	160690	255	1M623N1M	*	0	0	*	*
        """)

        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            self.make_ralign(Genome('gencode.v19'))

    @unittest.skipIf('CI' in os.environ, "avoid downloading full anno on CI")
    def test_detect_with_anno_resolve_with_format_library(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ENST00000466557.2.Intron6_neg	81	chr1	155831	255	1M8431N1M	*	0	0	*	*
        ENST00000496488.1.Intron1_pos	65	chr1	160690	255	1M623N1M	*	0	0	*	*
        """)

        gc19 = Genome('gencode.v19')
        self.make_ralign(gc19, library_format='ISF')
        with ReadAlignments(self.tmpout) as table:
            intron1 = gc19.interval('chr1', '+', 160690, 161313)
            reads1 = table.junctions.find_overlapping(intron1)
            self.assertEqual(len(reads1), 1)
            self.assertEqual(reads1[0].interval, intron1)

            intron6 = gc19.interval('chr1', '-', 155831, 164262)
            reads1 = table.junctions.find_overlapping(intron6)
            self.assertEqual(len(reads1), 1)
            self.assertEqual(reads1[0].interval, intron6)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full anno on CI")
    def test_detect_with_anno(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ENST00000456328.2.Intron1	0	chr1	12227	255	1M385N1M	*	0	0	*	*
        ENST00000438504.2.Intron1	0	chr1	24891	255	1M4429N1M	*	0	0	*	*
        """)

        gc19 = Genome('gencode.v19')
        positive_intron = gc19.interval("chr1", "+", 12227, 12612)
        negative_intron = gc19.interval("chr1", "-", 24891, 29320)

        self.make_ralign(gc19)
        with ReadAlignments(self.tmpout) as table:
            self.assertTrue(table.junctions.stranded)
            self.assertTrue(table.alignments.stranded)
            self.assertTrue(table.matches.stranded)
            self.assertFalse(table.variants.stranded)
            self.assertEqual(table.junctions[0].interval, positive_intron)
            self.assertEqual(table.junctions[1].interval, negative_intron)
            self.assertEqual(table.alignments[0].strand, '+')
            self.assertEqual(table.alignments[1].strand, '+')

        self.make_jralign(gc19)
        with JReadAlignments(self.tmpout) as table:
            self.assertTrue(table.junctions.stranded)
            self.assertEqual(table.junctions[0].interval, positive_intron)
            self.assertEqual(table.junctions[1].interval, negative_intron)
            self.assertEqual(table.junctions[0][0].strand, '+')
            self.assertEqual(table.junctions[1][0].strand, '+')

    def test_detect_with_library_format_unknown(self):
        dumptext(
            self.tmpsam, sam_header1, """
        Simple	0	chr1	0	255	1M1N1M	*	0	0	*	*
        """)

        with self.assertRaisesRegex(ValueError, 'library format'):
            self.make_ralign(library_format='unknown')

    def test_detect_with_library_format_unpaired(self):
        dumptext(self.tmpsam, sam_header1, """
        Unpaired	0	chr1	0	255	1M1N1M	*	0	0	*	*
        """)

        with self.assertRaisesRegex(ValueError, 'single end read'):
            self.make_ralign(library_format='ISF')

    def test_detect_with_library_format_stranded_forward(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ForwardFirst	65	chr1	0	255	1M1N1M	*	0	0	*	*
        ForwardSecond	129	chr1	10	255	1M1N1M	*	0	0	*	*
        ReverseFirst	81	chr1	20	255	1M1N1M	*	0	0	*	*
        ReverseSecond	145	chr1	30	255	1M1N1M	*	0	0	*	*
        """)
        self.make_ralign(library_format='ISF')
        with ReadAlignments(self.tmpout) as table:
            self.assertTrue(table.junctions.stranded)
            self.assertTrue(table.alignments.stranded)
            self.assertTrue(table.matches.stranded)
            self.assertFalse(table.variants.stranded)
            self.assertEqual(table.junctions[0].strand, '+')
            self.assertEqual(table.junctions[1].strand, '-')
            self.assertEqual(table.junctions[2].strand, '-')
            self.assertEqual(table.junctions[3].strand, '+')
            self.assertEqual(table.alignments[0].strand, '+')
            self.assertEqual(table.alignments[1].strand, '+')
            self.assertEqual(table.alignments[2].strand, '-')
            self.assertEqual(table.alignments[3].strand, '-')

        self.make_jralign(library_format='ISF')
        with JReadAlignments(self.tmpout) as table:
            self.assertTrue(table.junctions.stranded)
            self.assertEqual(table.junctions[0].strand, '+')
            self.assertEqual(table.junctions[1].strand, '-')
            self.assertEqual(table.junctions[2].strand, '-')
            self.assertEqual(table.junctions[3].strand, '+')
            self.assertEqual(table.junctions[0][0].strand, '+')
            self.assertEqual(table.junctions[1][0].strand, '+')
            self.assertEqual(table.junctions[2][0].strand, '-')
            self.assertEqual(table.junctions[3][0].strand, '-')

    def test_detect_with_library_format_stranded_reverse(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ForwardFirst	65	chr1	0	255	1M1N1M	*	0	0	*	*
        ForwardSecond	129	chr1	10	255	1M1N1M	*	0	0	*	*
        ReverseFirst	81	chr1	20	255	1M1N1M	*	0	0	*	*
        ReverseSecond	145	chr1	30	255	1M1N1M	*	0	0	*	*
        """)

        self.make_ralign(library_format='ISR')
        with ReadAlignments(self.tmpout) as table:
            self.assertEqual(table.junctions[0].strand, '-')
            self.assertEqual(table.junctions[1].strand, '+')
            self.assertEqual(table.junctions[2].strand, '+')
            self.assertEqual(table.junctions[3].strand, '-')
            self.assertEqual(table.alignments[0].strand, '+')
            self.assertEqual(table.alignments[1].strand, '+')
            self.assertEqual(table.alignments[2].strand, '-')
            self.assertEqual(table.alignments[3].strand, '-')

        self.make_jralign(library_format='ISR')
        with JReadAlignments(self.tmpout) as table:
            self.assertEqual(table.junctions[0].strand, '-')
            self.assertEqual(table.junctions[1].strand, '+')
            self.assertEqual(table.junctions[2].strand, '+')
            self.assertEqual(table.junctions[3].strand, '-')
            self.assertEqual(table.junctions[0][0].strand, '+')
            self.assertEqual(table.junctions[1][0].strand, '+')
            self.assertEqual(table.junctions[2][0].strand, '-')
            self.assertEqual(table.junctions[3][0].strand, '-')

    def test_rdist(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ForwardFirst	65	chr1	0	255	1M1N1M	*	0	0	*	*
        ForwardSecond	129	chr1	10	255	1M1N1M	*	0	0	*	*
        ReverseFirst	81	chr1	20	255	1M1N1M	*	0	0	*	*
        ReverseSecond	145	chr1	30	255	1M1N1M	*	0	0	*	*
        """)

        fd, self.tmprdist = tempfile.mkstemp()
        os.close(fd)

        self.make_jralign(library_format='ISF')
        ReadDistributions.build_rdist(self.tmprdist, [self.tmpout])
        with ReadDistributions(self.tmprdist) as table:
            self.assertTrue(table.junctions.stranded)
            self.assertEqual(table.junctions[0].strand, '+')
            self.assertEqual(table.junctions[1].strand, '-')
            self.assertEqual(table.junctions[2].strand, '-')
            self.assertEqual(table.junctions[3].strand, '+')

        os.unlink(self.tmprdist)

    @unittest.skipIf('CI' in os.environ, "avoid downloading full 2bit on CI")
    def test_pileup(self):
        dumptext(
            self.tmpsam, sam_header1, """
        ForwardFirstA	65	chr1	100001	255	1M	*	0	0	A	*
        ReverseSecondA	145	chr1	100001	255	1M	*	0	0	A	*
        ForwardFirstC	65	chr1	100001	255	1M	*	0	0	C	*
        """)
        genome = Genome('hg19')
        self.make_ralign(genome)
        with ReadAlignments(self.tmpout) as table:
            np.testing.assert_equal([[2, 1, 0, 0, 0]], table.pileup(Interval('chr1', '+', 100000, 100001, 'hg19')))

        self.make_ralign(genome, library_format='ISF')
        with ReadAlignments(self.tmpout) as table:
            np.testing.assert_equal([[1, 1, 0, 0, 0]], table.pileup(Interval('chr1', '+', 100000, 100001, 'hg19')))
            np.testing.assert_equal([[1, 0, 0, 0, 0]], table.pileup(Interval('chr1', '-', 100000, 100001, 'hg19')))


########################################################

if __name__ == "__main__":
    unittest.main()
