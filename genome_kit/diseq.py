from dataclasses import dataclass
from typing import Optional
from typing import Sequence, Literal

from genome_kit import Interval, Genome, Transcript

@dataclass
class _DisjointIntervalMetadata:
    transcript_id: str
    reference_genome: str
    chromosome: str
    transcript_strand: Literal["+", "-"]


class DisjointIntervalSequence:
    """
    A flattened 5′→3′ coordinate system over a sequence of disjoint exonic Intervals.

    The coordinate system is 0-based, half-open: [0, length).
    """


    def __init__(
            self,
            _intervals: Sequence[Interval],
            *,
            _metadata: _DisjointIntervalMetadata,
    ):
        """
        Low-level constructor.

        Parameters are not part of the stable public API. Use
        `from_transcript()` or `from_exons()` instead.
        """
        if len(_intervals) > 0:
            ival0 = _intervals[0]
            # check that the _intervals have the same chr and strand
            assert all([ival.chromosome == ival0.chromosome for ival in _intervals])
            assert all([ival.strand == ival0.strand for ival in _intervals])

            # check that the _intervals are consistently anchored
            assert all([ival.anchor == ival0.anchor for ival in _intervals])
            assert all([ival.anchor_offset == ival0.anchor_offset for ival in _intervals])
            # TODO reverse sort for negative strand?
            _intervals = sorted(_intervals, key=lambda ival: ival.start if ival.strand == "+" else -ival.end)
            self._intervals = _intervals
        else:
            self._intervals = []

        self.transcript_id = _metadata.transcript_id
        self.reference_genome = _metadata.reference_genome
        self.chromosome = _metadata.chromosome
        self.transcript_strand = _metadata.transcript_strand

    @classmethod
    def from_exons(
            cls,
            exons: Sequence[Interval],
            *,
            transcript_id: str | None = None,
            reference_genome: str | Genome | None = None,
            name: str | None = None,
    ) -> "DisjointIntervalSequence":

        return cls(exons, _metadata=_DisjointIntervalMetadata(
            transcript_id=transcript_id,
            reference_genome=reference_genome,
            chromosome=exons[0].chromosome,
            transcript_strand=exons[0].strand,
        ))

    @classmethod
    def from_transcript(
            cls,
            transcript: Transcript,
            *,
            region: Literal["exons", "cds", "utr5", "utr3"] = "exons",
            name: str | None = None,
    ) -> "DisjointIntervalSequence":
        if region == "exons":
            intervals = transcript.exons
        elif region == "cds":
            intervals = transcript.cds_exons
        elif region == "utr5":
            intervals = transcript.utr5_exons
        elif region == "utr3":
            intervals = transcript.utr3_exons
        else:
            raise ValueError(f"Invalid region: {region}")

        return cls(intervals, _metadata=_DisjointIntervalMetadata(
            transcript_id=transcript.id,
            reference_genome=transcript.reference_genome,
            chromosome=intervals[0].chromosome,
            transcript_strand=intervals[0].strand,
        ))


    @property
    def intervals(self) -> tuple[Interval, ...]:
        """Underlying genomic exonic intervals, sorted 5′→3′."""
        return tuple(self._intervals)

    @property
    def genomic_span(self) -> Interval:
        """Smallest Interval spanning all exons."""
        return Interval(
            chromosome=self._intervals[0].chromosome,
            start=self._intervals[0].start,
            end=self._intervals[-1].end,
            strand=self._intervals[0].strand,
            anchor=self._intervals[0].anchor,
            anchor_offset=self._intervals[0].anchor_offset,
        )

    @property
    def start(self) -> int:
        """Always 0 in DIS coordinates."""
        return 0

    @property
    def end(self) -> int:
        """Total length in bases of all exons combined."""
        return sum(i.length for i in self._intervals)

    @property
    def length(self) -> int:
        return self.end

    def __len__(self) -> int:
        return self.end

    @property
    def disjoint_start(self) -> int:
        """Always 0 in DIS coordinates."""
        return 0

    @property
    def disjoint_end(self) -> int:
        """Total length in bases of all exons combined."""
        return self.end


    # def shift(self, amount: int) -> "DisjointIntervalSequence":
    #     """
    #     Shift upstream/downstream in DIS coordinates.
    #
    #     Positive amount shifts towards the 3′ end of the transcript;
    #     negative shifts towards the 5′ end.
    #     """
    #
    # def expand(
    #         self,
    #         upstream: int,
    #         dnstream: int | None = None,
    # ) -> "DisjointIntervalSequence":
    #     """
    #     Expand in DIS coordinates. 'upstream' and 'dnstream' are 5′/3′
    #     relative to the transcript, regardless of genomic strand.
    #     """
    #
    def intersect(
            self,
            other: "DisjointIntervalSequence | Interval",
    ) -> "DisjointIntervalSequence | None":
        """
        Computes the intersection of the current interval sequence with another interval
        or disjoint interval sequence.

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval` | :py:class:`~genome_kit.DisjointIntervalSequence`
            The other interval or interval sequence to intersect with.

        Returns
        -------
        :py:class:`~genome_kit.DisjointIntervalSequence` | :data:`None`
            A DisjointIntervalSequence representing the intersection of intervals,
                 or `None` if there is no overlap.
        """
        if isinstance(other, Interval):
            intersected_intervals = []
            for interval in self.intervals:
                intersection = interval.intersect(other)
                if intersection is not None:
                    intersected_intervals.append(intersection)

            if not intersected_intervals:
                return None

            return DisjointIntervalSequence(intersected_intervals, _metadata=_DisjointIntervalMetadata(
                transcript_id=self.transcript_id,
                reference_genome=self.reference_genome,
                chromosome=self.chromosome,
                transcript_strand=self.transcript_strand,
            ))

        if isinstance(other, DisjointIntervalSequence):
            intersected_intervals = []
            for self_interval in self.intervals:
                for other_interval in other.intervals:
                    intersection = self_interval.intersect(other_interval)
                    if intersection is not None:
                        intersected_intervals.append(intersection)

            if not intersected_intervals:
                return None

            return DisjointIntervalSequence(intersected_intervals, _metadata=_DisjointIntervalMetadata(
                transcript_id=self.transcript_id,
                reference_genome=self.reference_genome,
                chromosome=self.chromosome,
                transcript_strand=self.transcript_strand,
            ))

        raise TypeError(f"Cannot intersect with type {type(other).__name__}")

    # def subtract(
    #         self,
    #         other: "DisjointIntervalSequence | Interval",
    # ) -> Sequence["DisjointIntervalSequence"]: ...
    #
    # # Predicates in DIS space (bool results)
    # def overlaps(self, other: "DisjointIntervalSequence | Interval") -> bool: ...
    # def contains(self, other: "DisjointIntervalSequence | Interval") -> bool: ...
    # def within(self, other: "DisjointIntervalSequence | Interval") -> bool: ...
    # def upstream_of(self, other: "DisjointIntervalSequence | Interval") -> bool: ...
    # def dnstream_of(self, other: "DisjointIntervalSequence | Interval") -> bool: ...
    #
    # def distance(
    #         self,
    #         other: "DisjointIntervalSequence | Interval | Sequence[Interval] | Sequence[DisjointIntervalSequence]",
    #         *,
    #         method: Literal["midpoint", "end5", "end3"] = "midpoint",
    # ) -> float | list[float]: ...


    def lift_from_transcript(
        self,
        src: Interval,
        *,
        clip: bool = False,
    ) -> "Optional[DisjointIntervalSequence]":
        """
        Lift a genomic/transcript Interval onto this DIS.

        Returns a DisjointIntervalSequence representing the interval in the DIS coordinate system.

        If clip=False:
            - raises a ValueError if `interval` does not lie fully within exons
              covered by this DIS.

        If clip=True:
            - returns the overlap of `interval` with this DIS projected
              into DIS coordinates, or None if there is no overlap at all.
        """
        if (src.chromosome != self.chromosome or src.strand != self.transcript_strand):
            raise ValueError(f"Interval chromosome/strand mismatch")

        # Find overlapping intervals and compute DIS coordinates
        result_intervals = []
        dis_offset = 0
        has_overlap = False

        for segment in self._intervals:
            overlap_start = max(src.start, segment.start)
            overlap_end = min(src.end, segment.end)

            if overlap_start < overlap_end:
                has_overlap = True

                if self.transcript_strand == "+":
                    exon_relative_start = overlap_start - segment.start
                    exon_relative_end = overlap_end - segment.start
                else:
                    exon_relative_start = segment.end - overlap_end
                    exon_relative_end = segment.end - overlap_start

                dis_start = dis_offset + exon_relative_start
                dis_end = dis_offset + exon_relative_end

                dis_interval = Interval(
                    reference_genome=self.reference_genome,
                    chromosome=self.chromosome,
                    strand=self.transcript_strand,
                    start=dis_start,
                    end=dis_end,
                    anchor=segment.anchor,
                    anchor_offset=segment.anchor_offset
                )
                result_intervals.append(dis_interval)

            dis_offset += len(segment)

        if not has_overlap:
            if clip:
                return None
            else:
                raise ValueError(f"Interval does not overlap with any exons")

        # Check if the input interval is fully contained (only when clip=False)
        if not clip:
            # Calculate total overlap length
            total_overlap = sum(len(interval) for interval in result_intervals)
            if total_overlap < len(src):
                raise ValueError(f"Interval extends beyond exon boundaries")

        # Create and return a new DisjointIntervalSequence with the result intervals
        return DisjointIntervalSequence(
            _intervals=result_intervals,
            _metadata=_DisjointIntervalMetadata(
                transcript_id=self.transcript_id,
                reference_genome=self.reference_genome,
                chromosome=self.chromosome,
                transcript_strand=self.transcript_strand,
            )
        )

    def lower_to_transcript(
        self,
        start: int,
        end: int,
    ) -> list[Interval]:
        """
        Lower [start, end) in DIS coordinates to one or more genomic
        Intervals in transcript/genomic coordinates.

        The returned Intervals:
            - lie on the same chromosome/strand/reference_genome as the exons,
            - are ordered 5′→3′ along the transcript,
            - are clipped to the specified start/end range.
        """
        if start < 0 or end < 0 or start >= end:
            return []

        result = []
        current_dis_pos = 0

        # Iterate through exons to find which ones overlap with [start, end)
        for segment in self._intervals:
            # Calculate the DIS coordinates for this exon
            seg_dis_start = current_dis_pos
            seg_dis_end = current_dis_pos + len(segment)

            # Check if this exon overlaps with the requested range
            if seg_dis_end > start and seg_dis_start < end:
                overlap_start = max(start, seg_dis_start)
                overlap_end = min(end, seg_dis_end)

                # Convert the overlap back to genomic coordinates within this exon
                exon_offset_start = overlap_start - seg_dis_start
                exon_offset_end = overlap_end - seg_dis_start

                # Create an interval for this portion of the exon
                # The exact implementation depends on how exons are represented
                # and the Interval class structure
                genomic_start = segment.start + exon_offset_start
                genomic_end = segment.start + exon_offset_end

                # Create interval with same chromosome, strand, and reference genome as exons
                interval = Interval(
                    chromosome=segment.chromosome,
                    start=genomic_start,
                    end=genomic_end,
                    strand=segment.strand,
                    reference_genome=segment.reference_genome
                )
                result.append(interval)

            current_dis_pos = seg_dis_end

            if current_dis_pos >= end:
                break

        return result
