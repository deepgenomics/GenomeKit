import unittest
from genome_kit import Interval, Genome
from genome_kit.diseq import DisjointIntervalSequence, _DisjointCoordinateMetadata


class TestLiftFromTranscript(unittest.TestCase):
    # TODO add negative strand tests for all cases

    def setUp(self):
        self.genome = Genome("hg19")

        exons_plus = [
            Interval("chr1", "+", 100, 200, self.genome.refg),
            Interval("chr1", "+", 300, 400, self.genome.refg),
        ]
        self.dis_plus = DisjointIntervalSequence(
            _intervals=exons_plus,
            _coord_metadata=_DisjointCoordinateMetadata(
                id="test_transcript_plus",
                reference_genome="hg19",
                chromosome="chr1",
                transcript_strand="+",
            )
        )

        exons_minus = [
            Interval("chr1", "-", 300, 400, self.genome.refg),  # First in 5'->3' order
            Interval("chr1", "-", 100, 200, self.genome.refg),  # Second in 5'->3' order
        ]
        self.dis_minus = DisjointIntervalSequence(
            _intervals=exons_minus,
            _coord_metadata=_DisjointCoordinateMetadata(
                id="test_transcript_minus",
                reference_genome="hg19",
                chromosome="chr1",
                transcript_strand="-",
            )
        )

    def test_chromosome_mismatch(self):
        interval = Interval("chr2", "+", 150, 180, self.genome.refg)

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=False)
        self.assertIn("chromosome/strand mismatch", str(cm.exception))

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=True)
        self.assertIn("chromosome/strand mismatch", str(cm.exception))

    def test_strand_mismatch(self):
        interval = Interval("chr1", "-", 150, 180, self.genome.refg)

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=False)
        self.assertIn("chromosome/strand mismatch", str(cm.exception))

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=True)
        self.assertIn("chromosome/strand mismatch", str(cm.exception))

    def test_no_overlap_no_clip(self):
        interval = Interval("chr1", "+", 250, 280, self.genome.refg)  # In intron

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=False)
        self.assertIn("does not overlap with any exons", str(cm.exception))

    def test_no_overlap_with_clip(self):
        interval = Interval("chr1", "+", 250, 280, self.genome.refg)  # In intron

        result = self.dis_plus.lift_from_transcript(interval, clip=True)
        self.assertIsNone(result)

    def test_fully_contained_single_exon_plus_strand(self):
        interval = Interval("chr1", "+", 120, 150, self.genome.refg)  # Within first exon

        result = self.dis_plus.lift_from_transcript(interval, clip=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 1)
        self.assertEqual(result.intervals[0].start, 20)  # 120 - 100 = 20 in DIS coordinates
        self.assertEqual(result.intervals[0].end, 50)    # 150 - 100 = 50 in DIS coordinates
        self.assertEqual(result.intervals[0].chromosome, "chr1")
        self.assertEqual(result.intervals[0].strand, "+")

    def test_fully_contained_single_exon_minus_strand(self):
        interval = Interval("chr1", "-", 320, 350, self.genome.refg)  # Within first exon

        result = self.dis_minus.lift_from_transcript(interval, clip=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 1)
        # For minus strand: exon is 300-400, interval is 320-350
        # In DIS coordinates: 400-350=50 to 400-320=80
        self.assertEqual(result.intervals[0].start, 50)
        self.assertEqual(result.intervals[0].end, 80)

    def test_spans_two_exons_plus_strand(self):
        interval = Interval("chr1", "+", 150, 350, self.genome.refg)  # Spans both exons

        result = self.dis_plus.lift_from_transcript(interval, clip=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 2)

        # First interval: 150-200 in genomic -> 50-100 in DIS
        self.assertEqual(result.intervals[0].start, 50)
        self.assertEqual(result.intervals[0].end, 100)

        # Second interval: 300-350 in genomic -> 100-150 in DIS
        self.assertEqual(result.intervals[1].start, 100)
        self.assertEqual(result.intervals[1].end, 150)

    def test_partial_overlap_with_clip(self):
        interval = Interval("chr1", "+", 50, 350, self.genome.refg)  # Entire first exon + part of second

        result = self.dis_plus.lift_from_transcript(interval, clip=True)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 2)

        # First interval: entire first exon 100-200 -> 0-100 in DIS
        self.assertEqual(result.intervals[0].start, 0)
        self.assertEqual(result.intervals[0].end, 100)

        # Second interval: partial second exon 300-350 -> 100-150 in DIS
        self.assertEqual(result.intervals[1].start, 100)
        self.assertEqual(result.intervals[1].end, 150)

    def test_partial_overlap_no_clip_raises(self):
        interval = Interval("chr1", "+", 50, 150, self.genome.refg)  # Extends before first exon

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=False)
        self.assertIn("extends beyond exon boundaries", str(cm.exception))

    def test_spans_intron_with_clip(self):
        interval = Interval("chr1", "+", 180, 320, self.genome.refg)  # Spans intron 200-300

        result = self.dis_plus.lift_from_transcript(interval, clip=True)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 2)

        # First part: 180-200 -> 80-100 in DIS
        self.assertEqual(result.intervals[0].start, 80)
        self.assertEqual(result.intervals[0].end, 100)

        # Second part: 300-320 -> 100-120 in DIS
        self.assertEqual(result.intervals[1].start, 100)
        self.assertEqual(result.intervals[1].end, 120)

    def test_spans_intron_no_clip_raises(self):
        interval = Interval("chr1", "+", 180, 320, self.genome.refg)  # Spans intron

        with self.assertRaises(ValueError) as cm:
            self.dis_plus.lift_from_transcript(interval, clip=False)
        self.assertIn("extends beyond exon boundaries", str(cm.exception))

    def test_exact_exon_boundaries(self):
        interval = Interval("chr1", "+", 100, 200, self.genome.refg)  # Exact first exon

        result = self.dis_plus.lift_from_transcript(interval, clip=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 1)
        self.assertEqual(result.intervals[0].start, 0)
        self.assertEqual(result.intervals[0].end, 100)

    def test_metadata_preservation(self):
        interval = Interval("chr1", "+", 120, 150, self.genome.refg)

        result = self.dis_plus.lift_from_transcript(interval, clip=False)

        self.assertEqual(result.id, "test_transcript_plus")
        self.assertEqual(result.reference_genome, "hg19")
        self.assertEqual(result.chromosome, "chr1")
        self.assertEqual(result.transcript_strand, "+")

    def test_empty_interval(self):
        interval = Interval("chr1", "+", 150, 150, self.genome.refg)

        result = self.dis_plus.lift_from_transcript(interval, clip=True)
        self.assertIsNone(result)

    def test_result_interval_properties(self):
        interval = Interval("chr1", "+", 120, 150, self.genome.refg)

        result = self.dis_plus.lift_from_transcript(interval, clip=False)
        result_interval = result.intervals[0]

        # Check that anchor properties are preserved from source exon
        source_exon = self.dis_plus.intervals[0]
        self.assertEqual(result_interval.anchor, source_exon.anchor)
        self.assertEqual(result_interval.anchor_offset, source_exon.anchor_offset)

    def test_single_base_interval(self):
        interval = Interval("chr1", "+", 120, 121, self.genome.refg)  # Single base

        result = self.dis_plus.lift_from_transcript(interval, clip=False)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 1)
        self.assertEqual(result.intervals[0].start, 20)
        self.assertEqual(result.intervals[0].end, 21)
        self.assertEqual(len(result.intervals[0]), 1)


class TestLowerToTranscript(unittest.TestCase):

    def setUp(self):
        # Create a simple DisjointIntervalSequence for testing
        # Create three exons on chromosome 1, positive strand
        # Exon 1: 100-200 (length 100, DIS 0-100)
        # Exon 2: 300-350 (length 50, DIS 100-150)
        # Exon 3: 500-600 (length 100, DIS 150-250)
        # Total DIS length: 250
        self.exons_positive = [
            Interval(
                chromosome="chr1",
                start=100,
                end=200,
                strand="+",
                reference_genome="hg38"
            ),
            Interval(
                chromosome="chr1",
                start=300,
                end=350,
                strand="+",
                reference_genome="hg38"
            ),
            Interval(
                chromosome="chr1",
                start=500,
                end=600,
                strand="+",
                reference_genome="hg38"
            )
        ]

        metadata_positive = _DisjointCoordinateMetadata(
            id="test_transcript",
            reference_genome="hg38",
            chromosome="chr1",
            transcript_strand="+",
        )

        self.dis_positive = DisjointIntervalSequence(self.exons_positive, _coord_metadata=metadata_positive)

        self.exons_negative = [
            Interval(
                chromosome="chr2",
                start=100,
                end=200,
                strand="-",
                reference_genome="hg38"
            ),
            Interval(
                chromosome="chr2",
                start=300,
                end=400,
                strand="-",
                reference_genome="hg38"
            )
        ]

        metadata_negative = _DisjointCoordinateMetadata(
            id="test_transcript_neg",
            reference_genome="hg38",
            chromosome="chr2",
            transcript_strand="-",
        )

        self.dis_negative = DisjointIntervalSequence(self.exons_negative, _coord_metadata=metadata_negative)

    def test_invalid_input_parameters(self):
        # Negative start
        result = self.dis_positive.lower_to_transcript(-1, 10)
        self.assertEqual(result, [])

        # Negative end
        result = self.dis_positive.lower_to_transcript(0, -1)
        self.assertEqual(result, [])

        # start >= end
        result = self.dis_positive.lower_to_transcript(10, 10)
        self.assertEqual(result, [])

        result = self.dis_positive.lower_to_transcript(10, 5)
        self.assertEqual(result, [])

    def test_single_exon_full_coverage(self):
        # Map DIS coordinates 10-20 (within first exon)
        result = self.dis_positive.lower_to_transcript(10, 20)

        self.assertEqual(len(result), 1)
        interval = result[0]

        # Should map to genomic coordinates 110-120
        self.assertEqual(interval.chromosome, "chr1")
        self.assertEqual(interval.start, 110)
        self.assertEqual(interval.end, 120)
        self.assertEqual(interval.strand, "+")
        self.assertEqual(interval.reference_genome, "hg38")

    def test_single_exon_partial_coverage(self):
        # Map DIS coordinates 50-80 (within first exon)
        result = self.dis_positive.lower_to_transcript(50, 80)

        self.assertEqual(len(result), 1)
        interval = result[0]

        # Should map to genomic coordinates 150-180
        self.assertEqual(interval.start, 150)
        self.assertEqual(interval.end, 180)

    def test_multiple_exons_coverage(self):
        # Map DIS coordinates 50-180 (spans all three exons)
        # Exon 1: DIS 0-100, genomic 100-200
        # Exon 2: DIS 100-150, genomic 300-350
        # Exon 3: DIS 150-250, genomic 500-600
        result = self.dis_positive.lower_to_transcript(50, 180)

        self.assertEqual(len(result), 3)

        # First interval (partial first exon: DIS 50-100 -> genomic 150-200)
        self.assertEqual(result[0].start, 150)
        self.assertEqual(result[0].end, 200)

        # Second interval (full second exon: DIS 100-150 -> genomic 300-350)
        self.assertEqual(result[1].start, 300)
        self.assertEqual(result[1].end, 350)

        # Third interval (partial third exon: DIS 150-180 -> genomic 500-530)
        self.assertEqual(result[2].start, 500)
        self.assertEqual(result[2].end, 530)

    def test_exon_boundary_mapping(self):
        # Map exactly the second exon (DIS 100-150)
        result = self.dis_positive.lower_to_transcript(100, 150)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 300)
        self.assertEqual(result[0].end, 350)

    def test_cross_exon_boundary(self):
        # Map DIS 90-110 (crosses first-second exon boundary)
        result = self.dis_positive.lower_to_transcript(90, 110)

        self.assertEqual(len(result), 2)

        # Part of first exon (DIS 90-100 -> genomic 190-200)
        self.assertEqual(result[0].start, 190)
        self.assertEqual(result[0].end, 200)

        # Part of second exon (DIS 100-110 -> genomic 300-310)
        self.assertEqual(result[1].start, 300)
        self.assertEqual(result[1].end, 310)

    def test_out_of_bounds_mapping(self):
        # Total DIS length is 250, try to map 200-300
        result = self.dis_positive.lower_to_transcript(200, 300)

        self.assertEqual(len(result), 1)
        # Should only cover the last exon from DIS 200-250 -> genomic 550-600
        self.assertEqual(result[0].start, 550)
        self.assertEqual(result[0].end, 600)

    def test_completely_out_of_bounds(self):
        # Total DIS length is 250, try to map 300-400
        result = self.dis_positive.lower_to_transcript(300, 400)

        # Should return empty list since range is completely out of bounds
        self.assertEqual(result, [])

    def test_zero_length_mapping(self):
        result = self.dis_positive.lower_to_transcript(0, 1)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 100)
        self.assertEqual(result[0].end, 101)

    def test_negative_strand_mapping(self):
        # Map DIS 10-20 on negative strand
        result = self.dis_negative.lower_to_transcript(10, 20)

        self.assertEqual(len(result), 1)
        interval = result[0]

        self.assertEqual(interval.chromosome, "chr2")
        self.assertEqual(interval.strand, "-")
        self.assertEqual(interval.reference_genome, "hg38")
        # For negative strand, should still map to genomic coordinates 110-120
        self.assertEqual(interval.start, 110)
        self.assertEqual(interval.end, 120)

    def test_result_ordering(self):
        # Map across multiple exons
        result = self.dis_positive.lower_to_transcript(50, 180)

        # Results should be in transcript order (5' to 3')
        self.assertEqual(len(result), 3)

        # For positive strand, genomic coordinates should be ascending
        self.assertLess(result[0].start, result[1].start)
        self.assertLess(result[1].start, result[2].start)
        # Non-overlapping intervals
        self.assertLessEqual(result[0].end, result[1].start)
        self.assertLessEqual(result[1].end, result[2].start)

    def test_metadata_preservation(self):
        result = self.dis_positive.lower_to_transcript(10, 20)

        self.assertEqual(len(result), 1)
        interval = result[0]

        self.assertEqual(interval.chromosome, self.dis_positive.chromosome)
        self.assertEqual(interval.strand, self.dis_positive.transcript_strand)
        self.assertEqual(interval.reference_genome, self.dis_positive.reference_genome)

    def test_edge_case_single_base(self):
        result = self.dis_positive.lower_to_transcript(0, 1)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 100)
        self.assertEqual(result[0].end, 101)
        self.assertEqual(len(result[0]), 1)

    def test_full_transcript_mapping(self):
        result = self.dis_positive.lower_to_transcript(0, 250)

        self.assertEqual(len(result), 3)

        # Should cover all three exons completely
        self.assertEqual(result[0].start, 100)
        self.assertEqual(result[0].end, 200)
        self.assertEqual(result[1].start, 300)
        self.assertEqual(result[1].end, 350)
        self.assertEqual(result[2].start, 500)
        self.assertEqual(result[2].end, 600)

    def test_gap_between_exons_mapping(self):
        # Map DIS 95-155 which spans across all three intervals but
        # includes gaps between them
        result = self.dis_positive.lower_to_transcript(95, 155)

        self.assertEqual(len(result), 3)

        # Should get parts of all three exons
        # First exon: DIS 95-100 -> genomic 195-200
        self.assertEqual(result[0].start, 195)
        self.assertEqual(result[0].end, 200)

        # Second exon: DIS 100-150 -> genomic 300-350 (full exon)
        self.assertEqual(result[1].start, 300)
        self.assertEqual(result[1].end, 350)

        # Third exon: DIS 150-155 -> genomic 500-505
        self.assertEqual(result[2].start, 500)
        self.assertEqual(result[2].end, 505)

    def test_start_at_exon_boundary(self):
        # Start at DIS 100 (beginning of second exon)
        result = self.dis_positive.lower_to_transcript(100, 120)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 300)
        self.assertEqual(result[0].end, 320)

    def test_end_at_exon_boundary(self):
        # End at DIS 100 (end of first exon)
        result = self.dis_positive.lower_to_transcript(80, 100)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 180)
        self.assertEqual(result[0].end, 200)


class TestDisjointIntervalSequenceIntersect(unittest.TestCase):

    def setUp(self):
        # Create basic intervals for testing
        self.interval1 = Interval('chr1', '+', 10, 20, 'hg19')
        self.interval2 = Interval('chr1', '+', 15, 25, 'hg19')
        self.interval3 = Interval('chr1', '+', 30, 40, 'hg19')
        self.interval4 = Interval('chr1', '+', 35, 45, 'hg19')

        # Non-overlapping intervals
        self.interval_far = Interval('chr1', '+', 100, 110, 'hg19')

        # Different chromosome
        self.interval_chr2 = Interval('chr2', '+', 10, 20, 'hg19')

        # Create DisjointIntervalSequences
        self.dis1 = DisjointIntervalSequence([self.interval1, self.interval3], _coord_metadata=_DisjointCoordinateMetadata(
            id="transcript1",
            reference_genome="hg19",
            chromosome="chr1",
            transcript_strand="+",
        ))
        self.dis2 = DisjointIntervalSequence([self.interval2, self.interval4], _coord_metadata=_DisjointCoordinateMetadata(
            id="transcript2",
            reference_genome="hg19",
            chromosome="chr1",
            transcript_strand="+",
        ))

    def test_intersect_with_disjoint_interval_sequence_full_overlap(self):
        # dis1: [10-20, 30-40]
        # dis2: [15-25, 35-45]
        result = self.dis1.intersect(self.dis2)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 2)
        self.assertEqual(result.intervals[0], Interval('chr1', '+', 15, 20, 'hg19'))
        self.assertEqual(result.intervals[1], Interval('chr1', '+', 35, 40, 'hg19'))

    def test_intersect_with_disjoint_interval_sequence_no_overlap(self):
        dis_no_overlap = DisjointIntervalSequence([self.interval_far], _coord_metadata=_DisjointCoordinateMetadata(
            id="transcript1",
            reference_genome="hg19",
            chromosome="chr1",
            transcript_strand="+",
        ))
        result = self.dis1.intersect(dis_no_overlap)

        self.assertIsNone(result)

    def test_intersect_with_disjoint_interval_sequence_partial_overlap(self):
        # Create sequence that overlaps only one interval in dis1
        partial_dis = DisjointIntervalSequence([Interval('chr1', '+', 5, 15, 'hg19')],
                                               _coord_metadata=_DisjointCoordinateMetadata(id="transcript1", reference_genome="hg19", chromosome="chr1", transcript_strand="+"))
        result = self.dis1.intersect(partial_dis)

        self.assertIsNotNone(result)
        self.assertEqual(len(result.intervals), 1)
        self.assertEqual(result.intervals[0], Interval('chr1', '+', 10, 15, 'hg19'))

    def test_intersect_empty_disjoint_sequence_with_interval(self):
        empty_dis = DisjointIntervalSequence([], _coord_metadata=_DisjointCoordinateMetadata(
            id="transcript1",
            reference_genome="hg19",
            chromosome="chr1",
            transcript_strand="+",
        ))
        result = empty_dis.intersect(self.interval1)

        self.assertIsNone(result)

    def test_intersect_commutative_with_same_metadata(self):
        result1 = self.dis1.intersect(self.dis2)
        result2 = self.dis2.intersect(self.dis1)

        # Both should have same intervals (though metadata may differ)
        if result1 is not None and result2 is not None:
            self.assertEqual(len(result1.intervals), len(result2.intervals))

    def test_intersect_empty_result_returns_none(self):
        non_overlapping_dis = DisjointIntervalSequence([Interval('chr1', '+', 50, 60, 'hg19')],
                                                       _coord_metadata=_DisjointCoordinateMetadata(id="transcript1", reference_genome="hg19", chromosome="chr1", transcript_strand="+"))
        result = self.dis1.intersect(non_overlapping_dis)

        self.assertIsNone(result)


class TestDisjointIntervalSequenceFromInterval(unittest.TestCase):

    def setUp(self):
        genome = Genome("gencode.v41")
        self.disjoint_intervals = [x.interval for x in genome.transcripts[2000].exons]
        self.transcript_id = genome.transcripts[2000].transcript_id
        self.other_transcript_id = genome.transcripts[2001].transcript_id

    def test_from_interval(self):
        dis = DisjointIntervalSequence.from_interval(self.disjoint_intervals, id=self.transcript_id)
        self.assertEqual(dis.coordinate_intervals, self.disjoint_intervals)
        self.assertEqual(dis.id, self.transcript_id)
        self.assertEqual(dis.reference_genome, self.disjoint_intervals[0].reference_genome)
        self.assertEqual(dis.chromosome, self.disjoint_intervals[0].chromosome)
        self.assertEqual(dis.transcript_strand, self.disjoint_intervals[0].strand)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, sum(len(interval) for interval in self.disjoint_intervals))

    def test_intervals_unordered(self):
        start, *middle, end = [x.interval for x in self.genome.transcripts[2000].exons]
        disordered_intervals = [end, *middle, start]
        disordered_dis = DisjointIntervalSequence.from_interval(disordered_intervals)
        ordered_dis = DisjointIntervalSequence.from_interval(self.disjoint_intervals)
        self.assertEqual(disordered_dis, ordered_dis)

    def test_no_transcript_id(self):
        dis = DisjointIntervalSequence.from_interval(self.disjoint_intervals)
        self.assertEqual(dis.coordinate_intervals, self.disjoint_intervals)
        self.assertEqual(dis.id, self.transcript_id)

    def test_mismatched_transcript_id(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_interval(
                self.disjoint_intervals, id=self.other_transcript_id)

    def test_nonexistent_transcript_id(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_interval(
                self.disjoint_intervals, id="invalid_transcript_id")

    def test_set_disjoint_interval_boundaries(self):
        dis = DisjointIntervalSequence.from_interval(
            self.disjoint_intervals, disjoint_interval_start=10, disjoint_interval_end=20)
        self.assertEqual(dis.start, 10)
        self.assertEqual(dis.end, 20)

    def test_set_disjoint_interval_boundaries_one_bp(self):
        dis = DisjointIntervalSequence.from_interval(
            self.disjoint_intervals, disjoint_interval_start=15, disjoint_interval_end=16)
        self.assertEqual(dis.start, 15)
        self.assertEqual(dis.end, 16)

    def test_disjoint_interval_boundaries_start_gt_end(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_interval(
                self.disjoint_intervals,
                disjoint_interval_start=20,
                disjoint_interval_end=10)

    def test_disjoint_interval_boundaries_negative_start(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_interval(
                self.disjoint_intervals,
                disjoint_interval_start=-5,
                disjoint_interval_end=10)

    def test_disjoint_interval_boundaries_zero_length(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_interval(
                self.disjoint_intervals,
                disjoint_interval_start=10,
                disjoint_interval_end=10)

    def test_disjoint_interval_boundaries_beyond_coord_end(self):
        coord_end = sum(len(interval) for interval in self.disjoint_intervals)
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_interval(
                self.disjoint_intervals,
                disjoint_interval_start=coord_end + 10,
                disjoint_interval_end=coord_end + 20)

    def test_empty_intervals(self):
        with self.assertRaises((ValueError)):
            DisjointIntervalSequence.from_interval([])

    def test_intervals_overlapping(self):
        ref_genome = self.disjoint_intervals[0].reference_genome
        chrom = self.disjoint_intervals[0].chromosome
        strand = self.disjoint_intervals[0].strand
        overlapping = [
            Interval(chrom, strand, 100, 200, ref_genome),
            Interval(chrom, strand, 150, 250, ref_genome),
        ]
        with self.assertRaises((ValueError)):
            DisjointIntervalSequence.from_interval(overlapping)

    def test_intervals_different_chromosomes(self):
        ref_genome = self.disjoint_intervals[0].reference_genome
        diff_chr = [
            Interval("chr1", "+", 100, 200, ref_genome),
            Interval("chr2", "+", 300, 400, ref_genome),
        ]
        with self.assertRaises((AssertionError, ValueError)):
            DisjointIntervalSequence.from_interval(diff_chr)

    def test_intervals_different_strands(self):
        ref_genome = self.disjoint_intervals[0].reference_genome
        chrom = self.disjoint_intervals[0].chromosome
        diff_strand = [
            Interval(chrom, "+", 100, 200, ref_genome),
            Interval(chrom, "-", 300, 400, ref_genome),
        ]
        with self.assertRaises((ValueError)):
            DisjointIntervalSequence.from_interval(diff_strand)

    def test_intervals_different_reference_genomes(self):
        diff_ref = [
            Interval("chr1", "+", 100, 200, "hg19"),
            Interval("chr1", "+", 300, 400, "hg38"),
        ]
        with self.assertRaises((ValueError)):
            DisjointIntervalSequence.from_interval(diff_ref)

    def test_single_interval(self):
        single = [self.disjoint_intervals[0]]
        dis = DisjointIntervalSequence.from_interval(single)
        self.assertEqual(len(dis.coordinate_intervals), 1)
        self.assertEqual(dis.coordinate_intervals[0], single[0])
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, len(single[0]))

    def test_adjacent_intervals(self):
        ref_genome = self.disjoint_intervals[0].reference_genome
        chrom = self.disjoint_intervals[0].chromosome
        strand = self.disjoint_intervals[0].strand
        adjacent_intervals = [
            Interval(chrom, strand, 100, 200, ref_genome),
            Interval(chrom, strand, 200, 300, ref_genome),
        ]
        dis = DisjointIntervalSequence.from_interval(adjacent_intervals)
        self.assertEqual(len(dis.coordinate_intervals), 2)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, 200)

    def test_only_disjoint_interval_start(self):
        dis = DisjointIntervalSequence.from_interval(
            self.disjoint_intervals, disjoint_interval_start=10)
        self.assertEqual(dis.start, 10)
        self.assertEqual(dis.end, sum(len(i) for i in self.disjoint_intervals))

    def test_only_disjoint_interval_end(self):
        dis = DisjointIntervalSequence.from_interval(
            self.disjoint_intervals, disjoint_interval_end=20)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, 20)


class TestDisjointIntervalSequenceFromTranscript(unittest.TestCase):

    def setUp(self):
        self.genome = Genome("gencode.v41")
        self.transcript = self.genome.transcripts[2000]
        self.disjoint_exons = [x.interval for x in self.transcript.exons]
        self.disjoint_cdss = [x.interval for x in self.transcript.cdss]
        self.disjoint_utr5s = [x.interval for x in self.transcript.utr5s]
        self.disjoint_utr3s = [x.interval for x in self.transcript.utr3s]
        self.transcript_id = self.transcript.transcript_id
        self.other_transcript_id = self.genome.transcripts[2001].transcript_id

    def test_from_transcript(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript)
        self.assertEqual(dis.coordinate_intervals, self.disjoint_exons)
        self.assertEqual(dis.id, self.transcript_id)
        self.assertEqual(dis.reference_genome, self.disjoint_exons[0].reference_genome)
        self.assertEqual(dis.chromosome, self.disjoint_exons[0].chromosome)
        self.assertEqual(dis.transcript_strand, self.disjoint_exons[0].strand)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, sum(len(interval) for interval in self.disjoint_exons))

    def test_from_transcript_region_exons(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="exons")
        self.assertEqual(dis.coordinate_intervals, self.disjoint_exons)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, sum(len(i) for i in self.disjoint_exons))

    def test_from_transcript_region_cds(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="cds")
        self.assertEqual(dis.coordinate_intervals, self.disjoint_cdss)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, sum(len(i) for i in self.disjoint_cdss))

    def test_from_transcript_region_utr5(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="utr5")
        self.assertEqual(dis.coordinate_intervals, self.disjoint_utr5s)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, sum(len(i) for i in self.disjoint_utr5s))

    def test_from_transcript_region_utr3(self):
        dis = DisjointIntervalSequence.from_transcript(self.transcript, region="utr3")
        self.assertEqual(dis.coordinate_intervals, self.disjoint_utr3s)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, sum(len(i) for i in self.disjoint_utr3s))

    def test_from_transcript_with_interval_boundaries(self):
        dis = DisjointIntervalSequence.from_transcript(
            self.transcript, disjoint_interval_start=10, disjoint_interval_end=20)
        self.assertEqual(dis.start, 10)
        self.assertEqual(dis.end, 20)
        self.assertEqual(dis.coordinate_intervals, self.disjoint_exons)

    def test_from_transcript_with_interval_boundaries_one_bp(self):
        dis = DisjointIntervalSequence.from_transcript(
            self.transcript, disjoint_interval_start=0, disjoint_interval_end=1)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, 1)
        self.assertEqual(dis.coordinate_intervals, self.disjoint_exons)

    def test_invalid_region_type(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(self.transcript, region="invalid")

    def test_invalid_region_type_introns(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(self.transcript, region="introns")

    def test_interval_boundaries_start_gt_end(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(
                self.transcript, disjoint_interval_start=20, disjoint_interval_end=10)

    def test_interval_boundaries_negative_start(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(
                self.transcript, disjoint_interval_start=-5, disjoint_interval_end=10)

    def test_interval_boundaries_zero_length(self):
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(
                self.transcript, disjoint_interval_start=10, disjoint_interval_end=10)

    def test_interval_boundaries_beyond_coord_end(self):
        coord_end = sum(len(i) for i in self.disjoint_exons)
        with self.assertRaises(ValueError):
            DisjointIntervalSequence.from_transcript(
                self.transcript,
                disjoint_interval_start=coord_end + 10,
                disjoint_interval_end=coord_end + 20)

    def test_regions_produce_distinct_intervals(self):
        dis_exons = DisjointIntervalSequence.from_transcript(self.transcript, region="exons")
        dis_cds = DisjointIntervalSequence.from_transcript(self.transcript, region="cds")
        dis_utr5 = DisjointIntervalSequence.from_transcript(self.transcript, region="utr5")
        dis_utr3 = DisjointIntervalSequence.from_transcript(self.transcript, region="utr3")

        # CDS, UTR5, UTR3 should each be distinct from full exons
        self.assertNotEqual(dis_exons.coordinate_intervals, dis_cds.coordinate_intervals)
        self.assertNotEqual(dis_exons.coordinate_intervals, dis_utr5.coordinate_intervals)
        self.assertNotEqual(dis_exons.coordinate_intervals, dis_utr3.coordinate_intervals)

        # UTRs and CDS should be mutually distinct
        self.assertNotEqual(dis_cds.coordinate_intervals, dis_utr5.coordinate_intervals)
        self.assertNotEqual(dis_cds.coordinate_intervals, dis_utr3.coordinate_intervals)
        self.assertNotEqual(dis_utr5.coordinate_intervals, dis_utr3.coordinate_intervals)

        # CDS + UTRs total length should equal exons total length
        self.assertEqual(dis_exons.end, dis_cds.end + dis_utr5.end + dis_utr3.end)

    def test_only_disjoint_interval_start(self):
        dis = DisjointIntervalSequence.from_transcript(
            self.transcript, disjoint_interval_start=10)
        self.assertEqual(dis.start, 10)
        self.assertEqual(dis.end, sum(len(i) for i in self.disjoint_exons))
        self.assertEqual(dis.coordinate_intervals, self.disjoint_exons)

    def test_only_disjoint_interval_end(self):
        dis = DisjointIntervalSequence.from_transcript(
            self.transcript, disjoint_interval_end=20)
        self.assertEqual(dis.start, 0)
        self.assertEqual(dis.end, 20)
        self.assertEqual(dis.coordinate_intervals, self.disjoint_exons)
