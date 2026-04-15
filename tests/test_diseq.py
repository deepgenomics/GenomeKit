import unittest
from genome_kit import Interval
from tests import MiniGenome
from genome_kit.diseq import DisjointIntervalSequence

REFG = "hg19.mini"


def _make_intervals(specs, refg = REFG):
    """Helper: specs is list of (chrom, strand, start, end)."""
    return [
        Interval(chrom, strand, start, end, refg) for chrom, strand, start, end in specs
    ]


class TestInit(unittest.TestCase):

    def test_non_interval_raises(self):
        with self.assertRaises(TypeError):
            DisjointIntervalSequence(["not an interval"])
        with self.assertRaises(TypeError):
            DisjointIntervalSequence([42])
        iv = Interval("chr1", "+", 100, 200, REFG)
        with self.assertRaises(TypeError):
            DisjointIntervalSequence([iv, "bad"])

    def test_empty_list_raises(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence([])

    def test_mixed_chromosomes_raises(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr2", "+", 300, 400)])
        with self.assertRaises(ValueError):
            DisjointIntervalSequence(ivs)

    def test_mixed_strands_raises(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "-", 300, 400)])
        with self.assertRaises(ValueError):
            DisjointIntervalSequence(ivs)

    def test_mixed_reference_genomes_raises(self):
        ivs = [
            Interval("chr2", "+", 100, 200, "hg38.p12.mini"),
            Interval("chr2", "+", 300, 400, "hg38.p13.mini"),
        ]
        with self.assertRaises(ValueError):
            DisjointIntervalSequence(ivs)

    def test_overlapping_intervals_raises(self):
        ivs = _make_intervals([("chr1", "+", 100, 250), ("chr1", "+", 200, 400)])
        with self.assertRaises(ValueError):
            DisjointIntervalSequence(ivs)

    def test_overlapping_intervals_negative_strand_raises(self):
        ivs = _make_intervals([("chr1", "-", 100, 250), ("chr1", "-", 200, 400)])
        with self.assertRaises(ValueError):
            DisjointIntervalSequence(ivs)

    def test_adjacent_intervals_ok(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 200, 300)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coordinate_length, 200)
        self.assertEqual(dis.coordinate_intervals, tuple(ivs))

    def test_sorts_out_of_order_positive(self):
        ivs = _make_intervals([("chr1", "+", 300, 400), ("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coordinate_intervals[0].start, 100)
        self.assertEqual(dis.coordinate_intervals[1].start, 300)

    def test_sorts_out_of_order_negative(self):
        ivs = _make_intervals([("chr1", "-", 100, 200), ("chr1", "-", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coordinate_intervals[0].end, 400)
        self.assertEqual(dis.coordinate_intervals[1].end, 200)

    def test_coordinate_length(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 450)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coordinate_length, 250)

    def test_start_end_default_to_full_interval(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, 200)

    def test_out_of_bounds_indices_allowed(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])  # coord_length=100
        dis = DisjointIntervalSequence(ivs, start=-10, end=50)
        self.assertEqual(dis.start, -10)
        dis2 = DisjointIntervalSequence(ivs, start=10, end=300)
        self.assertEqual(dis2.end, 300)
        dis3 = DisjointIntervalSequence(ivs, start=-100, end=300)
        self.assertEqual(dis3.start, -100)
        self.assertEqual(dis3.end, 300)
        dis4 = DisjointIntervalSequence(ivs, start=300, end=300)
        self.assertEqual(dis4.start, 300)
        self.assertEqual(dis4.end, 300)
        dis5 = DisjointIntervalSequence(ivs, start=-300, end=-300)
        self.assertEqual(dis5.start, -300)
        self.assertEqual(dis5.end, -300)

    def test_start_greater_than_end_raises(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        with self.assertRaises(ValueError):
            DisjointIntervalSequence(ivs, start=80, end=10)

    def test_start_equals_end_allowed(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, start=50, end=50)
        self.assertEqual(dis.length, 0)

    def test_custom_interval_indices(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, start=10, end=50)
        self.assertEqual(dis.start, 10)
        self.assertEqual(dis.end, 50)
        self.assertEqual(dis.length, 40)


class TestFromIntervals(unittest.TestCase):

    def setUp(self):
        self.genome = MiniGenome("gencode.v41")
        self.transcript = self.genome.transcripts["ENST00000233331.12"]

    def test_happy_path(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence.from_intervals(ivs, coord_name="mycoord")
        self.assertEqual(dis.coord_name, "mycoord")
        self.assertEqual(dis.coordinate_length, 200)
        self.assertEqual(dis.coordinate_intervals, tuple(ivs))

    def test_coord_and_interval_id_independent(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence.from_intervals(
            ivs, coord_name="c1", interval_name="i1"
        )
        self.assertEqual(dis.coord_name, "c1")
        self.assertEqual(dis.name, "i1")

    def test_single_interval(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence.from_intervals(ivs)
        self.assertEqual(len(dis.coordinate_intervals), 1)
        self.assertEqual(dis.coordinate_length, 100)
        self.assertEqual(dis.coordinate_intervals, tuple(ivs))

    def test_adjacent_intervals(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 200, 300)])
        dis = DisjointIntervalSequence.from_intervals(ivs)
        self.assertEqual(dis.coordinate_length, 200)
        self.assertEqual(dis.coordinate_intervals, tuple(ivs))

    def test_extracts_interval_from_exon(self):
        exons = list(self.transcript.exons)
        dis = DisjointIntervalSequence.from_intervals(exons)
        expected = tuple(e.interval for e in self.transcript.exons)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_extracts_interval_from_cds(self):
        cdss = list(self.transcript.cdss)
        dis = DisjointIntervalSequence.from_intervals(cdss)
        expected = tuple(c.interval for c in self.transcript.cdss)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_extracts_interval_from_utr5(self):
        utr5s = list(self.transcript.utr5s)
        dis = DisjointIntervalSequence.from_intervals(utr5s)
        expected = tuple(u.interval for u in self.transcript.utr5s)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_extracts_interval_from_utr3(self):
        utr3s = list(self.transcript.utr3s)
        dis = DisjointIntervalSequence.from_intervals(utr3s)
        expected = tuple(u.interval for u in self.transcript.utr3s)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_metadata_defaults_to_none(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence.from_intervals(ivs)
        self.assertEqual(dis.coord_name, None)
        self.assertEqual(dis.name, None)
        self.assertEqual(dis.reference_genome, REFG)
        self.assertEqual(dis.chromosome, "chr1")
        self.assertEqual(dis.coord_strand, "+")


class TestFromTranscript(unittest.TestCase):

    def setUp(self):
        self.genome = MiniGenome("gencode.v41")
        self.transcript = self.genome.transcripts["ENST00000233331.12"]
        self.neg_transcript = self.genome.transcripts["ENST00000448666.7"]

    def test_exons_region(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="exons")
        expected = tuple(e.interval for e in self.transcript.exons)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_cds_region(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="cds")
        expected = tuple(c.interval for c in self.transcript.cdss)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_utr5_region(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="utr5")
        expected = tuple(u.interval for u in self.transcript.utr5s)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_utr3_region(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="utr3")
        expected = tuple(u.interval for u in self.transcript.utr3s)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_metadata_defaults_to_transcript_id(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript)
        self.assertEqual(dis.coord_name, self.transcript.id)
        self.assertEqual(dis.name, self.transcript.id)
        self.assertEqual(dis.reference_genome, self.transcript.reference_genome)
        self.assertEqual(dis.chromosome, self.transcript.chromosome)
        self.assertEqual(dis.coord_strand, self.transcript.strand)

    def test_invalid_region_raises(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(self.transcript, region="invalid")

    def test_custom_id_overrides(self):
        dis = DisjointIntervalSequence.from_transcript(
            self.transcript, coord_name="custom_coord", interval_name="custom_iv"
        )
        self.assertEqual(dis.coord_name, "custom_coord")
        self.assertEqual(dis.name, "custom_iv")

    def test_negative_strand_exons(self):
        dis = DisjointIntervalSequence.from_transcript(self.neg_transcript, region="exons")
        expected = tuple(e.interval for e in self.neg_transcript.exons)
        self.assertEqual(dis.coordinate_intervals, expected)
        self.assertEqual(dis.coord_strand, "-")

    def test_negative_strand_cds(self):
        dis = DisjointIntervalSequence.from_transcript(self.neg_transcript, region="cds")
        expected = tuple(c.interval for c in self.neg_transcript.cdss)
        self.assertEqual(dis.coordinate_intervals, expected)

    def test_negative_strand_metadata(self):
        dis = DisjointIntervalSequence.from_transcript(self.neg_transcript)
        self.assertEqual(dis.coord_name, self.neg_transcript.id)
        self.assertEqual(dis.coord_strand, "-")
        self.assertEqual(dis.chromosome, self.neg_transcript.chromosome)


class TestProperties(unittest.TestCase):

    def test_metadata_getters_positive(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, coord_name="c", interval_name="i")
        self.assertEqual(dis.coord_name, "c")
        self.assertEqual(dis.name, "i")
        self.assertEqual(dis.reference_genome, REFG)
        self.assertEqual(dis.chromosome, "chr1")
        self.assertEqual(dis.coord_strand, "+")
        self.assertTrue(dis.on_coordinate_strand)

    def test_on_coordinate_strand_false(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertFalse(dis.on_coordinate_strand)

    def test_coord_transcript_strand_positive(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coord_strand, "+")

    def test_coord_transcript_strand_negative(self):
        ivs = _make_intervals([("chr1", "-", 100, 200)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coord_strand, "-")

    def test_coord_transcript_strand_unaffected_by_on_coordinate_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertEqual(dis.coord_strand, "+")

    def test_transcript_strand_on_coordinate_strand_positive(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertEqual(dis.strand, "+")

    def test_transcript_strand_on_coordinate_strand_negative(self):
        ivs = _make_intervals([("chr1", "-", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertEqual(dis.strand, "-")

    def test_transcript_strand_off_coordinate_strand_positive(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertEqual(dis.strand, "-")

    def test_transcript_strand_off_coordinate_strand_negative(self):
        ivs = _make_intervals([("chr1", "-", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertEqual(dis.strand, "+")

    def test_start_and_end(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, start=10, end=150)
        self.assertEqual(dis.start, 10)
        self.assertEqual(dis.end, 150)

    def test_coordinate_intervals(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, start=10, end=150)
        self.assertEqual(dis.coordinate_intervals, tuple(ivs))

    def test_coordinate_length(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.coordinate_length, 200)

    def test_length(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, start=10, end=150)
        self.assertEqual(dis.length, 140)

    def test_length_zero(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, start=50, end=50)
        self.assertEqual(dis.length, 0)


class TestStrandMethods(unittest.TestCase):

    def test_is_same_strand_true(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertTrue(dis.is_same_strand())

    def test_is_same_strand_false(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertFalse(dis.is_same_strand())

    def test_is_opposite_strand_true(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertTrue(dis.is_opposite_strand())

    def test_is_opposite_strand_false(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertFalse(dis.is_opposite_strand())

    def test_is_same_opposite_strand_mutually_exclusive(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis_same = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertNotEqual(dis_same.is_same_strand(), dis_same.is_opposite_strand())
        dis_opp = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertNotEqual(dis_opp.is_same_strand(), dis_opp.is_opposite_strand())

    def test_is_positive_strand_plus_on_coord(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertTrue(dis.is_positive_strand())

    def test_is_positive_strand_minus_off_coord(self):
        ivs = _make_intervals([("chr1", "-", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertTrue(dis.is_positive_strand())

    def test_is_positive_strand_false_plus_off_coord(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertFalse(dis.is_positive_strand())

    def test_is_positive_strand_false_minus_on_coord(self):
        ivs = _make_intervals([("chr1", "-", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertFalse(dis.is_positive_strand())

    def test_as_positive_strand_already_positive(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        result = dis.as_positive_strand()
        self.assertIs(result, dis)

    def test_as_positive_strand_flips(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(
            ivs, on_coordinate_strand=False, start=10, end=80
        )
        expected = DisjointIntervalSequence(
            ivs, on_coordinate_strand=True, start=10, end=80
        )
        result = dis.as_positive_strand()
        self.assertTrue(result.is_positive_strand())
        self.assertTrue(result.on_coordinate_strand)
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 80)
        self.assertEqual(result, expected)

    def test_as_negative_strand_already_negative(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        result = dis.as_negative_strand()
        self.assertIs(result, dis)

    def test_as_negative_strand_flips(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True, start=10, end=80)
        expected = DisjointIntervalSequence(
            ivs, on_coordinate_strand=False, start=10, end=80
        )
        result = dis.as_negative_strand()
        self.assertFalse(result.is_positive_strand())
        self.assertFalse(result.on_coordinate_strand)
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 80)
        self.assertEqual(result, expected)

    def test_flip_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        flipped = dis.flip_strand()
        self.assertFalse(flipped.on_coordinate_strand)
        flipped2 = flipped.flip_strand()
        self.assertTrue(flipped2.on_coordinate_strand)

    def test_flip_strand_preserves_coordinate_intervals(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        flipped = dis.flip_strand()
        self.assertEqual(flipped.coordinate_intervals, dis.coordinate_intervals)

    def test_flip_strand_preserves_start_end(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, start=10, end=80)
        flipped = dis.flip_strand()
        self.assertEqual(flipped.start, 10)
        self.assertEqual(flipped.end, 80)

    def test_flip_strand_preserves_metadata(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, coord_name="c", interval_name="i")
        flipped = dis.flip_strand()
        self.assertEqual(flipped.coord_name, "c")
        self.assertEqual(flipped.name, "i")

    def test_end5_end3_swap_on_flip_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, start=10, end=80)
        # On coordinate strand: end5 at start, end3 at end
        self.assertEqual(dis.end5_index, 10)
        self.assertEqual(dis.end3_index, 80)
        # Flipped: end5 at end, end3 at start
        flipped = dis.flip_strand()
        self.assertEqual(flipped.end5_index, 80)
        self.assertEqual(flipped.end3_index, 10)

    def test_as_opposite_strand_already_opposite(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        result = dis.as_opposite_strand()
        self.assertIs(result, dis)

    def test_as_opposite_strand_from_same(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True, start=10, end=80)
        result = dis.as_opposite_strand()
        self.assertFalse(result.on_coordinate_strand)
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 80)

    def test_as_opposite_strand_preserves_metadata(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, coord_name="c", interval_name="i")
        opp = dis.as_opposite_strand()
        self.assertEqual(opp.coord_name, "c")
        self.assertEqual(opp.name, "i")

    def test_as_same_strand_already_same(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        result = dis.as_same_strand()
        self.assertIs(result, dis)

    def test_as_same_strand_flips(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False, start=10, end=80)
        result = dis.as_same_strand()
        self.assertTrue(result.on_coordinate_strand)
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 80)

    def test_idempotency(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis_pos = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        self.assertIs(dis_pos.as_positive_strand().as_positive_strand(), dis_pos)
        self.assertIs(dis_pos.as_same_strand().as_same_strand(), dis_pos)
        dis_opp = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertIs(dis_opp.as_opposite_strand().as_opposite_strand(), dis_opp)

    def test_as_positive_strand_preserves_start_end(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False, start=10, end=80)
        result = dis.as_positive_strand()
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 80)

    def test_as_negative_strand_preserves_start_end(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=True, start=10, end=80)
        result = dis.as_negative_strand()
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 80)


class TestEndProperties(unittest.TestCase):

    def test_end5_default(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, coord_name="c", interval_name="i")
        # On coordinate strand: end5_index == start (0), end3_index == end (200)
        self.assertEqual(dis.end5_index, 0)
        self.assertEqual(dis.end3_index, 200)
        e5 = dis.end5
        self.assertEqual(len(e5), 0)
        self.assertEqual(e5.start, 0)
        self.assertEqual(e5.end, 0)
        self.assertEqual(e5.coord_name, "c")
        self.assertEqual(e5.name, None)
        expected = DisjointIntervalSequence(
            ivs, coord_name="c", interval_name=None, on_coordinate_strand=True, start=0, end=0
        )
        self.assertEqual(e5, expected)

    def test_end3_default(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        e3 = dis.end3
        self.assertEqual(len(e3), 0)
        self.assertEqual(e3.start, 200)
        self.assertEqual(e3.end, 200)
        self.assertEqual(dis.end3_index, 200)
        expected = DisjointIntervalSequence(
            ivs, on_coordinate_strand=True, start=200, end=200
        )
        self.assertEqual(e3, expected)

    def test_end5_end3_with_custom_indices(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, start=30, end=150)
        # On coordinate strand: end5_index == start, end3_index == end
        e5 = dis.end5
        e3 = dis.end3
        self.assertEqual(e5.start, 30)
        self.assertEqual(e5.end, 30)
        self.assertEqual(e3.start, 150)
        self.assertEqual(e3.end, 150)
        self.assertEqual(dis.end5_index, 30)
        self.assertEqual(dis.end3_index, 150)
        expected_e5 = DisjointIntervalSequence(
            ivs, on_coordinate_strand=True, start=30, end=30
        )
        expected_e3 = DisjointIntervalSequence(
            ivs, on_coordinate_strand=True, start=150, end=150
        )
        self.assertEqual(e5, expected_e5)
        self.assertEqual(e3, expected_e3)

    def test_coord_end5(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs)
        ce5 = dis.coord_end5
        self.assertEqual(len(ce5), 0)
        self.assertEqual(ce5.start, 0)
        self.assertEqual(ce5.end, 0)
        self.assertEqual(ce5.end5_index, 0)
        self.assertEqual(ce5.end3_index, 0)
        expected = DisjointIntervalSequence(
            ivs, on_coordinate_strand=True, start=0, end=0
        )
        self.assertEqual(ce5, expected)

    def test_coord_end3(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs)
        ce3 = dis.coord_end3
        self.assertEqual(len(ce3), 0)
        self.assertEqual(ce3.start, 100)
        self.assertEqual(ce3.end, 100)
        self.assertEqual(ce3.end5_index, 100)
        self.assertEqual(ce3.end3_index, 100)
        expected = DisjointIntervalSequence(
            ivs, on_coordinate_strand=True, start=100, end=100
        )
        self.assertEqual(ce3, expected)

    def test_end_preserves_coordinate_intervals(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(dis.end5.coordinate_intervals, dis.coordinate_intervals)
        self.assertEqual(dis.end3.coordinate_intervals, dis.coordinate_intervals)
        self.assertEqual(dis.coord_end5.coordinate_intervals, dis.coordinate_intervals)
        self.assertEqual(dis.coord_end3.coordinate_intervals, dis.coordinate_intervals)

    def test_end_preserves_on_coordinate_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertFalse(dis.end5.on_coordinate_strand)
        self.assertFalse(dis.end3.on_coordinate_strand)

    def test_coord_end_independent_of_on_coordinate_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertTrue(dis.coord_end5.on_coordinate_strand)
        self.assertTrue(dis.coord_end3.on_coordinate_strand)

    def test_end5_end3_opposite_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(
            ivs, start=30, end=150, on_coordinate_strand=False
        )
        # Off coord strand: end5 at end (150), end3 at start (30)
        self.assertEqual(dis.end5_index, 150)
        self.assertEqual(dis.end3_index, 30)
        e5 = dis.end5
        e3 = dis.end3
        self.assertEqual(e5.start, 150)
        self.assertEqual(e5.end, 150)
        self.assertEqual(e3.start, 30)
        self.assertEqual(e3.end, 30)


class TestDunderMethods(unittest.TestCase):

    def test_len(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs)
        self.assertEqual(len(dis), 200)

    def test_len_with_custom_indices(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs, start=10, end=80)
        self.assertEqual(len(dis), 70)

    def test_repr(self):
        ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
        dis = DisjointIntervalSequence(ivs, coord_name="ENST0001", interval_name="IV1")
        r = repr(dis)
        self.assertIn("DisjointIntervalSequence(", r)
        self.assertIn("coord_name='ENST0001'", r)
        self.assertIn("name='IV1'", r)
        self.assertIn("chr1:+", r)
        self.assertIn("len=200", r)
        self.assertIn('coord_intervals=(Interval("chr1", "+", 100, 200, "hg19.mini"), Interval("chr1", "+", 300, 400, "hg19.mini"))', r)
        self.assertIn("start=0", r)
        self.assertIn("end=200", r)
        self.assertIn("end5=0", r)
        self.assertIn("end3=200", r)

    def test_eq_same(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, coord_name="x", interval_name="i")
        b = DisjointIntervalSequence(ivs, coord_name="x", interval_name="i")
        self.assertEqual(a, b)

    def test_eq_different_coord_name(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, coord_name="x")
        b = DisjointIntervalSequence(ivs, coord_name="y")
        self.assertNotEqual(a, b)

    def test_eq_different_interval_name(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, interval_name="x")
        b = DisjointIntervalSequence(ivs, interval_name="y")
        self.assertNotEqual(a, b)

    def test_eq_different_on_coordinate_strand(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, on_coordinate_strand=True)
        b = DisjointIntervalSequence(ivs, on_coordinate_strand=False)
        self.assertNotEqual(a, b)

    def test_eq_different_chrom(self):
        a = DisjointIntervalSequence(_make_intervals([("chr1", "+", 100, 200)]))
        b = DisjointIntervalSequence(_make_intervals([("chr2", "+", 100, 200)]))
        self.assertNotEqual(a, b)

    def test_eq_different_refg(self):
        a = DisjointIntervalSequence([Interval("chr2", "+", 100, 200, "hg38.p12.mini")])
        b = DisjointIntervalSequence([Interval("chr2", "+", 100, 200, "hg38.p13.mini")])
        self.assertNotEqual(a, b)

    def test_eq_different_strand(self):
        a = DisjointIntervalSequence(_make_intervals([("chr1", "+", 100, 200)]))
        b = DisjointIntervalSequence(_make_intervals([("chr1", "-", 100, 200)]))
        self.assertNotEqual(a, b)

    def test_eq_different_coordinate_intervals(self):
        a = DisjointIntervalSequence(_make_intervals([("chr1", "+", 100, 200)]))
        b = DisjointIntervalSequence(_make_intervals([("chr1", "+", 100, 300)]))
        self.assertNotEqual(a, b)

    def test_eq_different_start(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, start=0, end=50)
        b = DisjointIntervalSequence(ivs, start=10, end=50)
        self.assertNotEqual(a, b)

    def test_eq_different_end(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, start=0, end=50)
        b = DisjointIntervalSequence(ivs, start=0, end=60)
        self.assertNotEqual(a, b)

    def test_eq_non_dis(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        dis = DisjointIntervalSequence(ivs)
        self.assertNotEqual(dis, "not a DIS")
        self.assertNotEqual(dis, 42)


# Helper for shift/expand/relational tests
# 2 exons on chr1+: [100,200) and [300,400), coordinate_length=200
_COORD_IVS = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 300, 400)])
_NEG_COORD_IVS = _make_intervals([("chr1", "-", 100, 200), ("chr1", "-", 300, 400)])


def _dis(
    start=0, end=200, on_coordinate_strand=True, coord_name="c", interval_name="i", ivs=None
):
    """Quick DIS factory for tests."""
    return DisjointIntervalSequence(
        ivs or _COORD_IVS,
        coord_name=coord_name,
        interval_name=interval_name,
        on_coordinate_strand=on_coordinate_strand,
        start=start,
        end=end,
    )


def _neg_dis(start=0, end=200, on_coordinate_strand=True):
    """Quick DIS factory for negative-strand coordinate interval tests."""
    return DisjointIntervalSequence(
        _NEG_COORD_IVS,
        coord_name="c",
        interval_name="i",
        on_coordinate_strand=on_coordinate_strand,
        start=start,
        end=end,
    )


class TestShift(unittest.TestCase):

    def test_shift_positive(self):
        dis = _dis(start=30, end=150)
        shifted = dis.shift(10)
        self.assertEqual(shifted.start, 40)
        self.assertEqual(shifted.end, 160)

    def test_shift_negative(self):
        dis = _dis(start=30, end=150)
        shifted = dis.shift(-10)
        self.assertEqual(shifted.start, 20)
        self.assertEqual(shifted.end, 140)

    def test_shift_zero(self):
        dis = _dis(start=30, end=150)
        shifted = dis.shift(0)
        self.assertEqual(shifted, dis)

    def test_shift_beyond_coordinate(self):
        dis = _dis(start=30, end=150)
        shifted = dis.shift(60)
        self.assertEqual(shifted.start, 90)
        self.assertEqual(shifted.end, 210)

    def test_shift_negative_beyond(self):
        dis = _dis(start=30, end=150)
        shifted = dis.shift(-40)
        self.assertEqual(shifted.start, -10)
        self.assertEqual(shifted.end, 110)

    def test_shift_zero_length(self):
        dis = _dis(start=50, end=50)
        shifted = dis.shift(5)
        self.assertEqual(shifted.start, 55)
        self.assertEqual(shifted.end, 55)
        self.assertEqual(shifted.length, 0)

    def test_shift_opposite_strand(self):
        # on_coordinate_strand=False: upstream_step=+1, downstream=-1
        # shift(10) downstream → subtract 10 from both
        dis = _dis(start=30, end=150, on_coordinate_strand=False)
        shifted = dis.shift(10)
        self.assertEqual(shifted.start, 20)
        self.assertEqual(shifted.end, 140)

    def test_shift_opposite_strand_negative_shift(self):
        # on_coordinate_strand=False: upstream_step=+1, downstream=-1
        # shift 10 upstream → add 10 to both
        dis = _dis(start=30, end=150, on_coordinate_strand=False)
        shifted = dis.shift(-10)
        self.assertEqual(shifted.start, 40)
        self.assertEqual(shifted.end, 160)

    def test_shift_preserves_metadata(self):
        dis = _dis(start=30, end=150, coord_name="mycoord", interval_name="myiv")
        shifted = dis.shift(10)
        self.assertEqual(shifted.coord_name, "mycoord")
        self.assertEqual(shifted.name, "myiv")
        self.assertTrue(shifted.on_coordinate_strand)

    def test_shift_preserves_metadata_opposite_strand(self):
        dis = _dis(
            start=30, end=150, coord_name="mycoord", interval_name="myiv",
            on_coordinate_strand=False,
        )
        shifted = dis.shift(10)
        self.assertEqual(shifted.coord_name, "mycoord")
        self.assertEqual(shifted.name, "myiv")
        self.assertFalse(shifted.on_coordinate_strand)

    def test_shift_preserves_coordinate_intervals(self):
        dis = _dis(start=30, end=150)
        shifted = dis.shift(10)
        self.assertEqual(shifted.coordinate_intervals, dis.coordinate_intervals)

    def test_shift_negative_strand_coords(self):
        dis = _neg_dis(start=30, end=150)
        shifted = dis.shift(10)
        self.assertEqual(shifted.start, 40)
        self.assertEqual(shifted.end, 160)

    def test_shift_negative_strand_coords_opposite(self):
        dis = _neg_dis(start=30, end=150, on_coordinate_strand=False)
        shifted = dis.shift(10)
        self.assertEqual(shifted.start, 20)
        self.assertEqual(shifted.end, 140)


class TestExpand(unittest.TestCase):

    def test_expand_symmetric(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(5)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 155)

    def test_expand_asymmetric(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(5, 10)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 160)

    def test_expand_upstream_only(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(5, 0)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 150)

    def test_expand_downstream_only(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(0, 10)
        self.assertEqual(expanded.start, 30)
        self.assertEqual(expanded.end, 160)

    def test_expand_zero(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(0)
        self.assertEqual(expanded, dis)

    def test_expand_negative_contracts(self):
        dis = _dis(start=30, end=150)
        contracted = dis.expand(-5, -10)
        self.assertEqual(contracted.start, 35)
        self.assertEqual(contracted.end, 140)

    def test_expand_contract_to_zero_length(self):
        dis = _dis(start=30, end=150)  # length=120
        contracted = dis.expand(-60, -60)
        self.assertEqual(contracted.start, 90)
        self.assertEqual(contracted.end, 90)
        self.assertEqual(contracted.length, 0)

    def test_expand_over_contraction_raises(self):
        dis = _dis(start=30, end=150)  # length=120
        with self.assertRaises(ValueError):
            dis.expand(-70, -70)

    def test_expand_opposite_strand(self):
        # on_coordinate_strand=False: upstream_step=+1
        # end5=150, end3=30. expand(5): end5 moves to 155, end3 moves to 25
        # start=min(155,25)=25, end=max(155,25)=155
        dis = _dis(start=30, end=150, on_coordinate_strand=False)
        expanded = dis.expand(5)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 155)

    def test_expand_opposite_strand_upstream_only(self):
        # on_coordinate_strand=False: upstream_step=+1
        # end5=150, end3=30. expand(5, 0): end5 moves to 155, end3 stays at 30
        dis = _dis(start=30, end=150, on_coordinate_strand=False)
        expanded = dis.expand(5, 0)
        self.assertEqual(expanded.start, 30)
        self.assertEqual(expanded.end, 155)

    def test_expand_opposite_strand_downstream_only(self):
        # on_coordinate_strand=False: upstream_step=+1
        # end5=150, end3=30. expand(0, 5): end5 stays at 150, end3 moves to 25
        dis = _dis(start=30, end=150, on_coordinate_strand=False)
        expanded = dis.expand(0, 5)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 150)

    def test_expand_zero_length_interval(self):
        dis = _dis(start=50, end=50)
        expanded = dis.expand(5)
        self.assertEqual(expanded.start, 45)
        self.assertEqual(expanded.end, 55)
        self.assertEqual(expanded.length, 10)

    def test_expand_beyond_coordinate(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(50, 0)
        self.assertEqual(expanded.start, -20)

    def test_expand_preserves_metadata(self):
        dis = _dis(start=30, end=150, coord_name="c", interval_name="i")
        expanded = dis.expand(5)
        self.assertEqual(expanded.coord_name, "c")
        self.assertEqual(expanded.name, "i")
        self.assertTrue(expanded.on_coordinate_strand)

    def test_expand_preserves_metadata_opposite_strand(self):
        dis = _dis(
            start=30, end=150, coord_name="c", interval_name="i",
            on_coordinate_strand=False,
        )
        expanded = dis.expand(5)
        self.assertEqual(expanded.coord_name, "c")
        self.assertEqual(expanded.name, "i")
        self.assertFalse(expanded.on_coordinate_strand)

    def test_expand_preserves_coordinate_intervals(self):
        dis = _dis(start=30, end=150)
        expanded = dis.expand(5)
        self.assertEqual(expanded.coordinate_intervals, dis.coordinate_intervals)

    def test_expand_negative_strand_coords(self):
        dis = _neg_dis(start=30, end=150)
        expanded = dis.expand(5, 10)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 160)

    def test_expand_negative_strand_coords_opposite(self):
        dis = _neg_dis(start=30, end=150, on_coordinate_strand=False)
        expanded = dis.expand(5)
        self.assertEqual(expanded.start, 25)
        self.assertEqual(expanded.end, 155)


class TestUpstreamOf(unittest.TestCase):

    def test_upstream_of_true(self):
        a = _dis(start=10, end=30)
        b = _dis(start=50, end=80)
        self.assertTrue(a.upstream_of(b))

    def test_upstream_of_false_overlap(self):
        a = _dis(start=10, end=60)
        b = _dis(start=50, end=80)
        self.assertFalse(a.upstream_of(b))

    def test_upstream_of_adjacent(self):
        a = _dis(start=10, end=50)
        b = _dis(start=50, end=80)
        self.assertTrue(a.upstream_of(b))

    def test_upstream_of_same_false(self):
        a = _dis(start=30, end=50)
        self.assertFalse(a.upstream_of(a))

    def test_upstream_of_zero_length(self):
        a = _dis(start=30, end=30)
        b = _dis(start=50, end=80)
        self.assertTrue(a.upstream_of(b))

    def test_upstream_of_zero_length_same_pos(self):
        a = _dis(start=50, end=50)
        b = _dis(start=50, end=80)
        self.assertTrue(a.upstream_of(b))

    def test_upstream_of_both_zero_length_same_pos(self):
        a = _dis(start=50, end=50)
        b = _dis(start=50, end=50)
        self.assertFalse(a.upstream_of(b))

    def test_upstream_of_opposite_strand(self):
        # on_coordinate_strand=False: upstream_step=+1, upstream = higher indices
        # a.start(100) >= b.end(80) → True
        a = _dis(start=100, end=150, on_coordinate_strand=False)
        b = _dis(start=50, end=80, on_coordinate_strand=False)
        self.assertTrue(a.upstream_of(b))

    def test_different_coord_space_raises(self):
        a = _dis(start=10, end=30)
        other_ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 500, 600)])
        b = _dis(start=10, end=30, ivs=other_ivs)
        with self.assertRaises(ValueError):
            a.upstream_of(b)

    def test_different_coord_name_allowed(self):
        a = _dis(start=10, end=30, coord_name="a")
        b = _dis(start=50, end=80, coord_name="b")
        self.assertTrue(a.upstream_of(b))

    def test_different_on_coord_strand_raises(self):
        a = _dis(start=10, end=30, on_coordinate_strand=True)
        b = _dis(start=50, end=80, on_coordinate_strand=False)
        with self.assertRaises(ValueError):
            a.upstream_of(b)

    def test_non_dis_raises(self):
        a = _dis(start=10, end=30)
        with self.assertRaises(TypeError):
            a.upstream_of("not a DIS")

    def test_upstream_of_negative_strand_coords(self):
        a = _neg_dis(start=10, end=30)
        b = _neg_dis(start=50, end=80)
        self.assertTrue(a.upstream_of(b))

    def test_upstream_of_negative_strand_coords_false(self):
        a = _neg_dis(start=50, end=80)
        b = _neg_dis(start=10, end=30)
        self.assertFalse(a.upstream_of(b))


class TestDnstreamOf(unittest.TestCase):

    def test_dnstream_of_true(self):
        a = _dis(start=50, end=80)
        b = _dis(start=10, end=30)
        self.assertTrue(a.dnstream_of(b))

    def test_dnstream_of_false(self):
        a = _dis(start=10, end=30)
        b = _dis(start=50, end=80)
        self.assertFalse(a.dnstream_of(b))

    def test_dnstream_of_adjacent(self):
        a = _dis(start=50, end=80)
        b = _dis(start=10, end=50)
        self.assertTrue(a.dnstream_of(b))

    def test_dnstream_of_same_false(self):
        a = _dis(start=30, end=50)
        self.assertFalse(a.dnstream_of(a))

    def test_dnstream_of_both_zero_length_same_pos(self):
        a = _dis(start=50, end=50)
        b = _dis(start=50, end=50)
        self.assertFalse(a.dnstream_of(b))

    def test_dnstream_of_opposite_strand(self):
        # on_coordinate_strand=False: upstream_step=+1
        # downstream = lower indices. a.end(80) <= b.start(100) → True
        a = _dis(start=50, end=80, on_coordinate_strand=False)
        b = _dis(start=100, end=150, on_coordinate_strand=False)
        self.assertTrue(a.dnstream_of(b))

    def test_different_coord_space_raises(self):
        a = _dis(start=50, end=80)
        other_ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 500, 600)])
        b = _dis(start=10, end=30, ivs=other_ivs)
        with self.assertRaises(ValueError):
            a.dnstream_of(b)

    def test_different_on_coord_strand_raises(self):
        a = _dis(start=50, end=80, on_coordinate_strand=True)
        b = _dis(start=30, end=60, on_coordinate_strand=False)
        with self.assertRaises(ValueError):
            a.dnstream_of(b)

    def test_dnstream_of_negative_strand_coords(self):
        a = _neg_dis(start=50, end=80)
        b = _neg_dis(start=10, end=30)
        self.assertTrue(a.dnstream_of(b))

    def test_dnstream_of_negative_strand_coords_false(self):
        a = _neg_dis(start=10, end=30)
        b = _neg_dis(start=50, end=80)
        self.assertFalse(a.dnstream_of(b))


class TestWithin(unittest.TestCase):

    def test_within_true(self):
        a = _dis(start=30, end=50)
        b = _dis(start=10, end=80)
        self.assertTrue(a.within(b))

    def test_within_false(self):
        a = _dis(start=10, end=80)
        b = _dis(start=30, end=50)
        self.assertFalse(a.within(b))

    def test_within_self(self):
        a = _dis(start=30, end=50)
        self.assertTrue(a.within(a))

    def test_within_zero_length(self):
        a = _dis(start=50, end=50)
        b = _dis(start=10, end=80)
        self.assertTrue(a.within(b))

    def test_within_at_boundary(self):
        a = _dis(start=10, end=80)
        b = _dis(start=10, end=80)
        self.assertTrue(a.within(b))

    def test_within_zero_length_at_boundary(self):
        a = _dis(start=10, end=10)
        b = _dis(start=10, end=80)
        c = _dis(start=80, end=80)
        self.assertTrue(a.within(b))
        self.assertTrue(c.within(b))

    def test_within_zero_length_outside(self):
        a = _dis(start=5, end=5)
        b = _dis(start=10, end=80)
        self.assertFalse(a.within(b))

    def test_within_opposite_strand(self):
        a = _dis(start=80, end=120, on_coordinate_strand=False)
        b = _dis(start=50, end=150, on_coordinate_strand=False)
        self.assertTrue(a.within(b))

    def test_different_coord_space_raises(self):
        a = _dis(start=30, end=50)
        other_ivs = _make_intervals([("chr1", "+", 100, 200), ("chr1", "+", 500, 600)])
        b = _dis(start=10, end=80, ivs=other_ivs)
        with self.assertRaises(ValueError):
            a.within(b)

    def test_different_on_coord_strand_raises(self):
        a = _dis(start=30, end=50, on_coordinate_strand=True)
        b = _dis(start=10, end=80, on_coordinate_strand=False)
        with self.assertRaises(ValueError):
            a.within(b)

    def test_non_dis_raises(self):
        a = _dis(start=30, end=50)
        with self.assertRaises(TypeError):
            a.within("not a DIS")

    def test_within_negative_strand_coords(self):
        a = _neg_dis(start=30, end=50)
        b = _neg_dis(start=10, end=80)
        self.assertTrue(a.within(b))

    def test_within_negative_strand_coords_false(self):
        a = _neg_dis(start=10, end=80)
        b = _neg_dis(start=30, end=50)
        self.assertFalse(a.within(b))


if __name__ == "__main__":
    unittest.main()
