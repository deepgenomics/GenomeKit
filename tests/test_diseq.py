import unittest
from genome_kit import Interval, Genome
from genome_kit.diseq import DisjointIntervalSequence, _DisjointIntervalMetadata


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
            _metadata=_DisjointIntervalMetadata(
                transcript_id="test_transcript_plus",
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
            _metadata=_DisjointIntervalMetadata(
                transcript_id="test_transcript_minus",
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

        self.assertEqual(result.transcript_id, "test_transcript_plus")
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
    """Tests for the lower_to_transcript method."""

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

        metadata_positive = _DisjointIntervalMetadata(
            transcript_id="test_transcript",
            reference_genome="hg38",
            chromosome="chr1",
            transcript_strand="+",
        )

        self.dis_positive = DisjointIntervalSequence(self.exons_positive, _metadata=metadata_positive)

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

        metadata_negative = _DisjointIntervalMetadata(
            transcript_id="test_transcript_neg",
            reference_genome="hg38",
            chromosome="chr2",
            transcript_strand="-",
        )

        self.dis_negative = DisjointIntervalSequence(self.exons_negative, _metadata=metadata_negative)

    def test_invalid_input_parameters(self):
        """Test that invalid input parameters return empty list."""
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
        """Test mapping within a single exon."""
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
        """Test mapping that covers part of an exon."""
        # Map DIS coordinates 50-80 (within first exon)
        result = self.dis_positive.lower_to_transcript(50, 80)

        self.assertEqual(len(result), 1)
        interval = result[0]

        # Should map to genomic coordinates 150-180
        self.assertEqual(interval.start, 150)
        self.assertEqual(interval.end, 180)

    def test_multiple_exons_coverage(self):
        """Test mapping that spans multiple exons."""
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
        """Test mapping exactly at exon boundaries."""
        # Map exactly the second exon (DIS 100-150)
        result = self.dis_positive.lower_to_transcript(100, 150)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 300)
        self.assertEqual(result[0].end, 350)

    def test_cross_exon_boundary(self):
        """Test mapping that crosses exon boundaries."""
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
        """Test mapping beyond the DIS length."""
        # Total DIS length is 250, try to map 200-300
        result = self.dis_positive.lower_to_transcript(200, 300)

        self.assertEqual(len(result), 1)
        # Should only cover the last exon from DIS 200-250 -> genomic 550-600
        self.assertEqual(result[0].start, 550)
        self.assertEqual(result[0].end, 600)

    def test_completely_out_of_bounds(self):
        """Test mapping completely beyond DIS bounds."""
        # Total DIS length is 250, try to map 300-400
        result = self.dis_positive.lower_to_transcript(300, 400)

        # Should return empty list since range is completely out of bounds
        self.assertEqual(result, [])

    def test_zero_length_mapping(self):
        """Test mapping at the very start."""
        result = self.dis_positive.lower_to_transcript(0, 1)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 100)
        self.assertEqual(result[0].end, 101)

    def test_negative_strand_mapping(self):
        # self.exons_negative = [
        #     Interval( chromosome="chr2", start=100, end=200, strand="-", reference_genome="hg38" ),
        #     Interval( chromosome="chr2", start=300, end=400, strand="-", reference_genome="hg38" )
        # ]

        """Test mapping on negative strand."""
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
        """Test that results are ordered 5' to 3' along transcript."""
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
        """Test that chromosome, strand, and reference genome are preserved."""
        result = self.dis_positive.lower_to_transcript(10, 20)

        self.assertEqual(len(result), 1)
        interval = result[0]

        self.assertEqual(interval.chromosome, self.dis_positive.chromosome)
        self.assertEqual(interval.strand, self.dis_positive.transcript_strand)
        self.assertEqual(interval.reference_genome, self.dis_positive.reference_genome)

    def test_edge_case_single_base(self):
        """Test mapping a single base."""
        result = self.dis_positive.lower_to_transcript(0, 1)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 100)
        self.assertEqual(result[0].end, 101)
        self.assertEqual(len(result[0]), 1)

    def test_full_transcript_mapping(self):
        """Test mapping the entire transcript."""
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
        """Test mapping across gaps between exons."""
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
        """Test mapping starting exactly at an exon boundary."""
        # Start at DIS 100 (beginning of second exon)
        result = self.dis_positive.lower_to_transcript(100, 120)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 300)
        self.assertEqual(result[0].end, 320)

    def test_end_at_exon_boundary(self):
        """Test mapping ending exactly at an exon boundary."""
        # End at DIS 100 (end of first exon)
        result = self.dis_positive.lower_to_transcript(80, 100)

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].start, 180)
        self.assertEqual(result[0].end, 200)
