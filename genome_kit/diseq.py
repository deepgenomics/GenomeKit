from copy import deepcopy
import enum
from dataclasses import dataclass
from typing import Sequence, Literal

from .interval import Interval
from .genome_annotation import Transcript


class IndexDirection(enum.Enum):
    """Controls how indices are assigned to the 5' and 3' ends of a coordinate space.
    Changing this value will result in undefined behaviour for
    :py:class:`DisjointIntervalSequence` objects created under a different convention.

    ``TRANSCRIPT_FIVE_TO_THREE``
        Index 0 is always at the coordinate transcript's 5' end, regardless of genomic strand.

    ``POSITIVE_STRAND_LEFT_TO_RIGHT``
        Index 0 is always at the leftmost genomic position relative to the positive strand.
        On the negative strand, this means the 3' end is at index 0.
    """

    TRANSCRIPT_FIVE_TO_THREE = "transcript_five_to_three"
    POSITIVE_STRAND_LEFT_TO_RIGHT = "positive_strand_left_to_right"


@dataclass(frozen=True)
class _CoordinateMetadata:
    name: str | None
    reference_genome: str
    chromosome: str
    transcript_strand: Literal["+", "-"]


@dataclass(frozen=True)
class _SegmentMetadata:
    name: str | None
    on_coordinate_strand: bool


class DisjointIntervalSequence:
    """A flattened coordinate system over a sequence of disjoint genomic Intervals.

    A DIS represents two layers:

    - A **coordinate space** defined by a sequence of non-overlapping genomic
      :py:class:`~genome_kit.Interval` objects (e.g. the exons of a transcript),
      which are flattened into a contiguous 0-based index space. Indices for the
      coordinate space are assigned according to the current :py:class:`IndexDirection`
      value.
    - A **segment** within that coordinate space, defined by a 5' and 3' index.
      The segment may lie on the same, or opposite, strand as the coordinate space.

    Use :py:meth:`from_transcript` or :py:meth:`from_intervals` to construct
    instances rather than calling the constructor directly.
    """

    _index_direction: IndexDirection = IndexDirection.TRANSCRIPT_FIVE_TO_THREE

    @classmethod
    def set_index_direction(cls, direction: IndexDirection) -> None:
        """Set the index direction convention for all DIS instances.

        Parameters
        ----------
        direction : :py:class:`IndexDirection`
            The index direction convention to use.
        """
        cls._index_direction = direction

    @classmethod
    def get_index_direction(cls) -> IndexDirection:
        """Get the current index direction convention.

        Returns
        -------
        :py:class:`IndexDirection`
        """
        return cls._index_direction

    def __init__(
        self,
        coordinate_intervals: Sequence[Interval],
        *,
        coord_name: str | None = None,
        segment_name: str | None = None,
        on_coordinate_strand: bool = True,
        start: int | None = None,
        end: int | None = None,
    ):
        """Low-level constructor.

        Prefer :py:meth:`from_transcript` or :py:meth:`from_intervals` for
        public construction.

        Parameters
        ----------
        coordinate_intervals : Sequence[:py:class:`~genome_kit.Interval`]
            Non-empty sequence of non-overlapping Intervals on the same
            chromosome, strand, and reference genome.
        coord_name : :py:class:`str` or None
            Optional name for the coordinate space.
        segment_name : :py:class:`str` or None
            Optional name for the segment.
        on_coordinate_strand : :py:class:`bool`
            Whether the segment is on the same strand as the coordinate
            intervals. Defaults to True. Can be used to represent a sequence that binds
            to the transcript if set to False.
        start : :py:class:`int` or None
            start index of the segment in the coordinate space. Defaults to 0
        end : :py:class:`int` or None
            end index of the segment in the coordinate space. Defaults to the length
            of the coordinate space.

        Raises
        ------
        ValueError
            If coordinate intervals are empty, inconsistent, overlapping, or if start
            is greater than end.
        TypeError
            If any element is not an Interval.
        """
        if len(coordinate_intervals) == 0:
            raise ValueError("coordinate_intervals must be non-empty")

        for i, iv in enumerate(coordinate_intervals):
            if not isinstance(iv, Interval):
                raise TypeError(
                    f"coordinate_intervals[{i}] is {type(iv).__name__}, expected Interval"
                )
            if iv.anchor is not None:
                raise ValueError(
                    f"coordinate_intervals[{i}] has an anchor set; "
                    f"anchored Intervals are not supported"
                )

        # Consistent chromosome, strand, reference_genome
        iv0 = coordinate_intervals[0]
        for iv in coordinate_intervals[1:]:
            if iv.chromosome != iv0.chromosome:
                raise ValueError(
                    f"All intervals must share the same chromosome, "
                    f"got {iv0.chromosome!r} and {iv.chromosome!r}"
                )
            if iv.strand != iv0.strand:
                raise ValueError(
                    f"All intervals must share the same strand, "
                    f"got {iv0.strand!r} and {iv.strand!r}"
                )
            if iv.reference_genome != iv0.reference_genome:
                raise ValueError(
                    f"All intervals must share the same reference genome, "
                    f"got {iv0.reference_genome!r} and {iv.reference_genome!r}"
                )

        # Sort 5'->3'
        if iv0.strand == "+":
            sorted_intervals = sorted(coordinate_intervals, key=lambda iv: iv.start)
        else:
            # On negative strand end is the 5' end since start < end.
            # Sort by -end to get 5'->3' order.
            sorted_intervals = sorted(coordinate_intervals, key=lambda iv: -iv.end)

        # No overlaps (adjacent/touching OK)
        for i in range(len(sorted_intervals) - 1):
            cur_iv, next_iv = sorted_intervals[i], sorted_intervals[i + 1]
            plus_strand_overlap = iv0.strand == "+" and cur_iv.end > next_iv.start
            minus_strand_overlap = iv0.strand == "-" and cur_iv.start < next_iv.end
            if plus_strand_overlap or minus_strand_overlap:
                raise ValueError(
                    f"Intervals must not overlap: [{cur_iv.start}, {cur_iv.end}) and [{next_iv.start}, {next_iv.end})"
                )

        self._coordinate_intervals: tuple[Interval, ...] = tuple(sorted_intervals)

        self._coord_metadata = _CoordinateMetadata(
            name=coord_name,
            reference_genome=iv0.reference_genome,
            chromosome=iv0.chromosome,
            transcript_strand=iv0.strand,
        )
        self._segment_metadata = _SegmentMetadata(
            name=segment_name,
            on_coordinate_strand=on_coordinate_strand,
        )

        # Default segment start/end to span the full coordinate
        if start is None:
            start = 0
        if end is None:
            end = self.coordinate_length

        # Validate that start is less than or equal to end
        if start > end:
            raise ValueError(
                f"start index {start} cannot be greater than end index {end}"
            )

        self._start: int = start
        self._end: int = end

    @classmethod
    def from_intervals(
        cls,
        intervals: Sequence[Interval],
        *,
        coord_name: str | None = None,
        segment_name: str | None = None,
    ) -> "DisjointIntervalSequence":
        """Construct a DIS from a sequence of Intervals
        (or :py:class:`~genome_kit.Exon`/:py:class:`~genome_kit.Cds`/:py:class:`~genome_kit.Utr` objects).

        Parameters
        ----------
        intervals : Sequence[:py:class:`~genome_kit.Interval`]
            Sequence of Interval or annotation objects.
        coord_name : :py:class:`str` or None
            Optional name for the coordinate space.
        segment_name : :py:class:`str` or None
            Optional name for the segment.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        coord_intervals = [
            iv.interval if hasattr(iv, "interval") else iv for iv in intervals
        ]
        return cls(coord_intervals, coord_name=coord_name, segment_name=segment_name)

    @classmethod
    def from_transcript(
        cls,
        transcript: Transcript,
        *,
        region: Literal["exons", "cds", "utr5", "utr3"] = "exons",
        coord_name: str | None = None,
        segment_name: str | None = None,
    ) -> "DisjointIntervalSequence":
        """Construct a DIS from a transcript's exons, CDS, or UTR regions.

        Parameters
        ----------
        transcript : :py:class:`~genome_kit.Transcript`
            The source Transcript object.
        region : :py:class:`str`
            Which region to extract — ``"exons"``, ``"cds"``,
            ``"utr5"``, or ``"utr3"``.
        coord_name : :py:class:`str` or None
            Optional name for the coordinate space. Defaults to ``transcript.id``.
        segment_name : :py:class:`str` or None
            Optional name for the segment. Defaults to ``transcript.id``.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`

        Raises
        ------
        ValueError
            If region is not one of the allowed values.
        """
        match region:
            case "exons":
                region_elements = transcript.exons
            case "cds":
                region_elements = transcript.cdss
            case "utr5":
                region_elements = transcript.utr5s
            case "utr3":
                region_elements = transcript.utr3s
            case _:
                raise ValueError(f"Invalid region: {region!r}")
        coord_intervals = [element.interval for element in region_elements]
        if coord_name is None:
            coord_name = transcript.id
        if segment_name is None:
            segment_name = transcript.id

        return cls(coord_intervals, coord_name=coord_name, segment_name=segment_name)

    @property
    def coord_name(self) -> str | None:
        """Name of the coordinate space, or None."""
        return self._coord_metadata.name

    @property
    def reference_genome(self) -> str:
        """Reference genome of the coordinate intervals."""
        return self._coord_metadata.reference_genome

    @property
    def chromosome(self) -> str:
        """Chromosome of the coordinate intervals."""
        return self._coord_metadata.chromosome

    @property
    def coord_strand(self) -> Literal["+", "-"]:
        """Strand of the coordinate intervals (the transcript strand)."""
        return self._coord_metadata.transcript_strand

    @property
    def name(self) -> str | None:
        """Name of the segment, or None."""
        return self._segment_metadata.name

    @property
    def on_coordinate_strand(self) -> bool:
        """True if the segment is on the same strand as the coordinate intervals."""
        return self._segment_metadata.on_coordinate_strand

    @property
    def strand(self) -> Literal["+", "-"]:
        """Effective strand of the segment, accounting for on_coordinate_strand."""
        if self.on_coordinate_strand:
            return self.coord_strand
        # Segment is on opposite strand
        if self.coord_strand == "+":
            return "-"
        return "+"

    @property
    def coordinate_end5_index(self) -> int:
        """5' index of the coordinate space."""
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return 0
        if self.coord_strand == "+":
            return 0
        return self.coordinate_length

    @property
    def coordinate_end3_index(self) -> int:
        """3' index of the coordinate space."""
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return self.coordinate_length
        if self.coord_strand == "+":
            return self.coordinate_length
        return 0

    @property
    def end5_index(self) -> int:
        """5' index of the segment."""
        if self._upstream_index_step() == -1:
            return self._start
        return self._end

    @property
    def end3_index(self) -> int:
        """3' index of the segment."""
        if self._upstream_index_step() == -1:
            return self._end
        return self._start

    @property
    def start(self) -> int:
        """Start index of the segment in the coordinate space."""
        return self._start

    @property
    def end(self) -> int:
        """End index of the segment in the coordinate space."""
        return self._end

    def _at_index(
        self, idx: int, on_coordinate_strand: bool
    ) -> "DisjointIntervalSequence":
        """Return a 0-length DIS at the given index position."""
        return DisjointIntervalSequence(
            self._coordinate_intervals,
            coord_name=self._coord_metadata.name,
            on_coordinate_strand=on_coordinate_strand,
            start=idx,
            end=idx,
        )

    @property
    def end5(self) -> "DisjointIntervalSequence":
        """0-length DIS at the segment's 5' end."""
        return self._at_index(
            self.end5_index, on_coordinate_strand=self.on_coordinate_strand
        )

    @property
    def end3(self) -> "DisjointIntervalSequence":
        """0-length DIS at the segment's 3' end."""
        return self._at_index(
            self.end3_index, on_coordinate_strand=self.on_coordinate_strand
        )

    @property
    def coord_end5(self) -> "DisjointIntervalSequence":
        """0-length DIS at the coordinate space's 5' end."""
        return self._at_index(self.coordinate_end5_index, on_coordinate_strand=True)

    @property
    def coord_end3(self) -> "DisjointIntervalSequence":
        """0-length DIS at the coordinate space's 3' end."""
        return self._at_index(self.coordinate_end3_index, on_coordinate_strand=True)

    @property
    def coordinate_intervals(self) -> tuple[Interval, ...]:
        """The underlying genomic intervals of the coordinate-space, sorted 5'->3'."""
        # Deepcopy to preserve imutability of this DIS
        return deepcopy(self._coordinate_intervals)

    @property
    def coordinate_length(self) -> int:
        """Total length of the coordinate space in bases."""
        return sum(len(iv) for iv in self._coordinate_intervals)

    @property
    def length(self) -> int:
        """Length of the segment on the coordinate space."""
        return self.end - self.start

    def _upstream_index_step(self, on_coordinate_strand: bool | None = None) -> int:
        """Return +1 or -1 indicating the upstream direction in index space.

        Args:
            on_coordinate_strand: Override for which strand to compute the step for.
                Defaults to this segment's on_coordinate_strand.
        """
        if on_coordinate_strand is None:
            on_coordinate_strand = self.on_coordinate_strand
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return -1 if on_coordinate_strand else 1
        # POSITIVE_STRAND_LEFT_TO_RIGHT: effective strand determines direction
        return -1 if self.strand == "+" else 1

    def _validate_same_coordinate_space(
        self, other: "DisjointIntervalSequence"
    ) -> None:
        """Raise if other does not share the same coordinate space."""
        if not isinstance(other, DisjointIntervalSequence):
            raise TypeError(
                f"Expected DisjointIntervalSequence, got {type(other).__name__}"
            )
        if self._coordinate_intervals != other._coordinate_intervals:
            raise ValueError("DIS objects must share the same coordinate intervals")

    def _from_end_indices(self, end5: int, end3: int) -> "DisjointIntervalSequence":
        """Return a new DIS with the same coordinate space but different segment indices."""
        # Validate end5 is upstream of or equal to end3
        if self._upstream_index_step() == -1:
            if end5 > end3:
                raise ValueError(
                    f"Invalid indices: end5 index {end5} is downstream of end3 index {end3}"
                )
        if self._upstream_index_step() == 1:
            if end5 < end3:
                raise ValueError(
                    f"Invalid indices: end5 index {end5} is downstream of end3 index {end3}"
                )
        return DisjointIntervalSequence(
            self._coordinate_intervals,
            coord_name=self._coord_metadata.name,
            segment_name=self._segment_metadata.name,
            on_coordinate_strand=self.on_coordinate_strand,
            start=min(end5, end3),
            end=max(end5, end3),
        )

    def shift(self, amount: int) -> "DisjointIntervalSequence":
        """Shift the segment downstream by amount (negative shifts upstream).

        The coordinate space is unchanged. Only the segment indices move.
        """
        downstream_step = -self._upstream_index_step()
        delta = amount * downstream_step
        return self._from_end_indices(
            self.end5_index + delta,
            self.end3_index + delta,
        )

    def expand(
        self, upstream: int, dnstream: int | None = None
    ) -> "DisjointIntervalSequence":
        """Expand the segment upstream and/or downstream.

        Negative values contract the segment. Raises ValueError if contraction
        would result in end5 being downstream of end3.

        Args:
            upstream: Bases to expand (or contract if negative) toward the 5' end.
            dnstream: Bases to expand (or contract if negative) toward the 3' end.
                Defaults to upstream (symmetric).
        """
        if dnstream is None:
            dnstream = upstream
        up_step = self._upstream_index_step()
        down_step = -up_step
        new_end5 = self.end5_index + (upstream * up_step)
        new_end3 = self.end3_index + (dnstream * down_step)
        # Validate end5 is still upstream of or equal to end3
        if (new_end5 - new_end3) * up_step < 0:
            raise ValueError(
                "Invalid expansion: end5 would be downstream of end3 "
                f"(end5={new_end5}, end3={new_end3})"
            )
        return self._from_end_indices(new_end5, new_end3)

    def upstream_of(self, other: "DisjointIntervalSequence") -> bool:
        """True if self is strictly upstream of other (no overlap).

        Requires the same coordinate space and same on_coordinate_strand.
        """
        self._validate_same_coordinate_space(other)
        if self.on_coordinate_strand != other.on_coordinate_strand:
            raise ValueError(
                f"Cannot compare: self is on "
                f"{'same' if self.on_coordinate_strand else 'opposite'} "
                f"strand but other is on "
                f"{'same' if other.on_coordinate_strand else 'opposite'} strand"
            )
        if self.length == 0 and other.length == 0 and self.start == other.start:
            return False
        if self._upstream_index_step() == -1:
            return self._end <= other.start
        return self._start >= other.end

    def dnstream_of(self, other: "DisjointIntervalSequence") -> bool:
        """True if self is strictly downstream of other (no overlap).

        Requires the same coordinate space and same on_coordinate_strand.
        """
        self._validate_same_coordinate_space(other)
        if self.on_coordinate_strand != other.on_coordinate_strand:
            raise ValueError(
                f"Cannot compare: self is on "
                f"{'same' if self.on_coordinate_strand else 'opposite'} "
                f"strand but other is on "
                f"{'same' if other.on_coordinate_strand else 'opposite'} strand"
            )
        if self.length == 0 and other.length == 0 and self.start == other.start:
            return False
        if self._upstream_index_step() == -1:
            return self._start >= other.end
        return self._end <= other.start

    def within(self, other: "DisjointIntervalSequence") -> bool:
        """True if self's segment is contained within other's segment.

        Requires the same coordinate space and same on_coordinate_strand.
        """
        self._validate_same_coordinate_space(other)
        if self.on_coordinate_strand != other.on_coordinate_strand:
            raise ValueError(
                f"Cannot compare: self is on "
                f"{'same' if self.on_coordinate_strand else 'opposite'} "
                f"strand but other is on "
                f"{'same' if other.on_coordinate_strand else 'opposite'} strand"
            )
        return self._start >= other.start and self._end <= other.end

    def is_same_strand(self) -> bool:
        """True if the segment is on the same strand as the coordinate intervals.

        Returns
        -------
        :py:class:`bool`
        """
        return self.on_coordinate_strand

    def is_positive_strand(self) -> bool:
        """True if the segment is on the positive strand.

        Returns
        -------
        :py:class:`bool`
        """
        if self.strand == "+":
            return True
        return False

    def as_positive_strand(self) -> "DisjointIntervalSequence":
        """Return a DIS with the segment on the positive strand.

        Returns ``self`` if already on the positive strand. The coordinate
        intervals are unchanged; only the segment strand is affected.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        if self.is_positive_strand():
            return self
        return self.flip_strand()

    def as_negative_strand(self) -> "DisjointIntervalSequence":
        """Return a DIS with the segment on the negative strand.

        Returns ``self`` if already on the negative strand. The coordinate
        intervals are unchanged; only the segment strand is affected.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        if not self.is_positive_strand():
            return self
        return self.flip_strand()

    def as_opposite_strand(self) -> "DisjointIntervalSequence":
        """Return a DIS with the segment on the opposite strand.

        Returns ``self`` if already on the opposite strand. The coordinate
        intervals are unchanged; only the segment strand is affected.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        if not self.on_coordinate_strand:
            return self
        return self.flip_strand()

    def as_same_strand(self) -> "DisjointIntervalSequence":
        """Return a DIS with the segment on the coordinate strand.

        Returns ``self`` if already on the coordinate strand. The coordinate
        intervals are unchanged; only the segment strand is affected.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        if self.on_coordinate_strand:
            return self
        return self.flip_strand()

    def flip_strand(self) -> "DisjointIntervalSequence":
        """Return a new DIS with ``on_coordinate_strand`` toggled.

        The coordinate intervals are unchanged. The segment's
        ``on_coordinate_strand`` is flipped.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        return DisjointIntervalSequence(
            self._coordinate_intervals,
            coord_name=self._coord_metadata.name,
            segment_name=self._segment_metadata.name,
            on_coordinate_strand=not self.on_coordinate_strand,
            start=self._start,
            end=self._end,
        )

    def __len__(self) -> int:
        """Return the length of the segment."""
        return self.length

    def __repr__(self) -> str:
        """Return a human-readable representation."""
        return (
            f"DisjointIntervalSequence("
            f"coord_name={self._coord_metadata.name!r}, "
            f"name={self._segment_metadata.name!r}, "
            f"{self.chromosome}:{self.coord_strand}, "
            f"len={self.length}, "
            f"coord_intervals={self._coordinate_intervals}, "
            f"start={self._start}, "
            f"end={self._end}, "
            f"end5={self.end5_index}, "
            f"end3={self.end3_index})"
        )

    def __eq__(self, other: object) -> bool:
        """Equality based on coordinate intervals, metadata, and index values."""
        if not isinstance(other, DisjointIntervalSequence):
            return NotImplemented
        return (
            self._coord_metadata == other._coord_metadata
            and self._segment_metadata == other._segment_metadata
            and self._start == other._start
            and self._end == other._end
            and self._coordinate_intervals == other._coordinate_intervals
        )
