import unittest
from genome_kit import Interval, Genome
from genome_kit.diseq import DisjointIntervalSequence

REFG = "hg19"


def _make_intervals(specs, refg=REFG):
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
            Interval("chr1", "+", 100, 200, "hg19"),
            Interval("chr1", "+", 300, 400, "hg38"),
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
        self.genome = Genome("gencode.v41")
        self.transcript = self.genome.transcripts[2002]

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
        self.genome = Genome("gencode.v41")
        self.transcript = self.genome.transcripts[2002]

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
        self.assertIn('coord_intervals=(Interval("chr1", "+", 100, 200, "hg19"), Interval("chr1", "+", 300, 400, "hg19"))', r)
        self.assertIn("start=0", r)
        self.assertIn("end=200", r)
        self.assertIn("end5=0", r)
        self.assertIn("end3=200", r)

    def test_eq_same(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, coord_name="x", interval_name="i")
        b = DisjointIntervalSequence(ivs, coord_name="x", interval_name="i")
        self.assertEqual(a, b)

    def test_eq_different_coord_id(self):
        ivs = _make_intervals([("chr1", "+", 100, 200)])
        a = DisjointIntervalSequence(ivs, coord_name="x")
        b = DisjointIntervalSequence(ivs, coord_name="y")
        self.assertNotEqual(a, b)

    def test_eq_different_interval_id(self):
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
        a = DisjointIntervalSequence([Interval("chr1", "+", 100, 200, "hg19")])
        b = DisjointIntervalSequence([Interval("chr1", "+", 100, 200, "hg38")])
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


if __name__ == "__main__":
    unittest.main()
