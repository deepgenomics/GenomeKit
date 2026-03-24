import enum
from dataclasses import dataclass
from typing import Sequence, Literal

from genome_kit import Interval, Transcript


# ---------------------------------------------------------------------------
# Index direction toggle
# ---------------------------------------------------------------------------

class IndexDirection(enum.Enum):
    TRANSCRIPT_FIVE_TO_THREE = "transcript_five_to_three"
    POSITIVE_STRAND_LEFT_TO_RIGHT = "positive_strand_left_to_right"


# ---------------------------------------------------------------------------
# Metadata dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _CoordinateMetadata:
    id: str | None
    reference_genome: str
    chromosome: str
    transcript_strand: Literal["+", "-"]


@dataclass(frozen=True)
class _IntervalMetadata:
    id: str | None
    on_coordinate_strand: bool


# ---------------------------------------------------------------------------
# DisjointIntervalSequence
# ---------------------------------------------------------------------------

class DisjointIntervalSequence:
    """A flattened coordinate system over a sequence of disjoint genomic Intervals."""

    _index_direction: IndexDirection = IndexDirection.TRANSCRIPT_FIVE_TO_THREE

    @classmethod
    def set_index_direction(cls, direction: IndexDirection) -> None:
        """Set the index direction convention for all DIS instances."""
        cls._index_direction = direction

    @classmethod
    def get_index_direction(cls) -> IndexDirection:
        """Get the current index direction convention."""
        return cls._index_direction

    def __init__(
        self,
        coordinate_intervals: Sequence[Interval],
        *,
        coord_id: str | None = None,
        interval_id: str | None = None,
        on_coordinate_strand: bool = True,
        interval_end5_index: int | None = None,
        interval_end3_index: int | None = None,
    ):
        """
        Low-level constructor.

        Parameters are not part of the stable public API. Use
        `from_transcript()` or `from_exons()` instead.

        Args:
            coordinate_intervals: Non-empty sequence of non-overlapping Intervals
                on the same chromosome, strand, and reference genome.
            coord_id: Optional identifier for the coordinate space.
            interval_id: Optional identifier for the interval.
            on_coordinate_strand: Whether the interval is on the same strand as
                the coordinate intervals.
            interval_end5_index: 5' index of the interval in the coordinate space.
                Defaults to the full coordinate span.
            interval_end3_index: 3' index of the interval in the coordinate space.
                Defaults to the full coordinate span.

        Raises:
            ValueError: If intervals are empty, inconsistent, or overlapping,
                or if index values are out of range.
            TypeError: If any element is not an Interval.
        """
        if len(coordinate_intervals) == 0:
            raise ValueError("coordinate_intervals must be non-empty")

        for i, iv in enumerate(coordinate_intervals):
            if not isinstance(iv, Interval):
                raise TypeError(
                    f"coordinate_intervals[{i}] is {type(iv).__name__}, expected Interval"
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

        self._coordinate_length: int = sum(len(iv) for iv in sorted_intervals)

        self._coord_metadata = _CoordinateMetadata(
            id=coord_id,
            reference_genome=iv0.reference_genome,
            chromosome=iv0.chromosome,
            transcript_strand=iv0.strand,
        )
        self._interval_metadata = _IntervalMetadata(
            id=interval_id,
            on_coordinate_strand=on_coordinate_strand,
        )

        # Default interval end5/end3 to span the full coordinate
        if interval_end5_index is None:
            interval_end5_index = self.coordinate_end5_index if on_coordinate_strand else self.coordinate_end3_index
        if interval_end3_index is None:
            interval_end3_index = self.coordinate_end3_index if on_coordinate_strand else self.coordinate_end5_index

        # Validate that end5 is upstream of or equal to end3
        up_step = self._upstream_index_step(on_coordinate_strand)
        if (interval_end5_index - interval_end3_index) * up_step < 0:
            raise ValueError(
                f"interval_end5_index={interval_end5_index} must be upstream of or equal to "
                f"interval_end3_index={interval_end3_index} (upstream_step={up_step})"
            )

        self._interval_end5_index: int = interval_end5_index
        self._interval_end3_index: int = interval_end3_index

    # -------------------------------------------------------------------
    # Public constructors
    # -------------------------------------------------------------------

    @classmethod
    def from_intervals(
        cls,
        intervals: Sequence[Interval],
        *,
        coord_id: str | None = None,
        interval_id: str | None = None,
    ) -> "DisjointIntervalSequence":
        """Construct a DIS from a sequence of Intervals (or Exon/Cds/Utr objects).

        If elements have an ``.interval`` attribute (e.g. Exon, Cds, Utr),
        the plain Interval is extracted automatically.

        Args:
            intervals: Sequence of Interval or annotation objects with ``.interval``.
            coord_id: Optional identifier for the coordinate space.
            interval_id: Optional identifier for the interval.
        """
        # Extract .interval if items are Exon/Cds/Utr
        coord_intervals = []
        for iv in intervals:
            if type(iv) is not Interval and hasattr(iv, "interval"):
                coord_intervals.append(iv.interval)
            else:
                coord_intervals.append(iv)
        return cls(coord_intervals, coord_id=coord_id, interval_id=interval_id)

    @classmethod
    def from_transcript(
        cls,
        transcript: Transcript,
        *,
        region: Literal["exons", "cds", "utr5", "utr3"] = "exons",
        coord_id: str | None = None,
        interval_id: str | None = None,
    ) -> "DisjointIntervalSequence":
        """Construct a DIS from a transcript's exons, CDS, or UTR regions.

        Args:
            transcript: The source Transcript object.
            region: Which region to extract — "exons", "cds", "utr5", or "utr3".
            coord_id: Optional coordinate ID. Defaults to ``transcript.id``.
            interval_id: Optional interval ID. Defaults to ``transcript.id``.

        Raises:
            ValueError: If region is not one of the allowed values.
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
        if coord_id is None:
            coord_id = transcript.id
        if interval_id is None:
            interval_id = transcript.id

        return cls(coord_intervals, coord_id=coord_id, interval_id=interval_id)

    # -------------------------------------------------------------------
    # Coordinate metadata properties
    # -------------------------------------------------------------------

    @property
    def coord_id(self) -> str | None:
        """Identifier for the coordinate space, or None."""
        return self._coord_metadata.id

    @property
    def reference_genome(self) -> str:
        """Reference genome of the coordinate intervals."""
        return self._coord_metadata.reference_genome

    @property
    def chromosome(self) -> str:
        """Chromosome of the coordinate intervals."""
        return self._coord_metadata.chromosome

    @property
    def coord_transcript_strand(self) -> Literal["+", "-"]:
        """Strand of the coordinate intervals (the transcript strand)."""
        return self._coord_metadata.transcript_strand

    # -------------------------------------------------------------------
    # Interval metadata properties
    # -------------------------------------------------------------------

    @property
    def id(self) -> str | None:
        """Identifier for the interval, or None."""
        return self._interval_metadata.id

    @property
    def on_coordinate_strand(self) -> bool:
        """True if the interval is on the same strand as the coordinate intervals."""
        return self._interval_metadata.on_coordinate_strand

    @property
    def transcript_strand(self) -> Literal["+", "-"]:
        """Effective strand of the interval, accounting for on_coordinate_strand."""
        if self.on_coordinate_strand:
            return self.coord_transcript_strand
        # Interval is on opposite strand
        if self.coord_transcript_strand == "+":
            return "-"
        return "+"

    # -------------------------------------------------------------------
    # Coordinate index properties (toggle-dependent)
    # -------------------------------------------------------------------

    @property
    def coordinate_end5_index(self) -> int:
        """5' index of the coordinate space. Depends on the index direction toggle."""
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return 0
        # direction is POSITIVE_STRAND_LEFT_TO_RIGHT
        if self.coord_transcript_strand == "+":
            return 0
        return self._coordinate_length

    @property
    def coordinate_end3_index(self) -> int:
        """3' index of the coordinate space. Depends on the index direction toggle."""
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return self._coordinate_length
        if self.coord_transcript_strand == "+":
            return self._coordinate_length
        return 0

    # -------------------------------------------------------------------
    # Interval index properties
    # -------------------------------------------------------------------

    @property
    def interval_end5_index(self) -> int:
        """5' index of the interval in the coordinate space."""
        return self._interval_end5_index

    @property
    def interval_end3_index(self) -> int:
        """3' index of the interval in the coordinate space."""
        return self._interval_end3_index

    # -------------------------------------------------------------------
    # End properties (0-length DIS at boundaries)
    # -------------------------------------------------------------------

    def _at_index(self, idx: int, on_coordinate_strand: bool) -> "DisjointIntervalSequence":
        """Return a 0-length DIS at the given index position."""
        return DisjointIntervalSequence(
            self._coordinate_intervals,
            coord_id=self._coord_metadata.id,
            on_coordinate_strand=on_coordinate_strand,
            interval_end5_index=idx,
            interval_end3_index=idx,
        )

    @property
    def end5(self) -> "DisjointIntervalSequence":
        """0-length DIS at the interval's 5' end."""
        return self._at_index(self._interval_end5_index, on_coordinate_strand=self.on_coordinate_strand)

    @property
    def end3(self) -> "DisjointIntervalSequence":
        """0-length DIS at the interval's 3' end."""
        return self._at_index(self._interval_end3_index, on_coordinate_strand=self.on_coordinate_strand)

    @property
    def coord_end5(self) -> "DisjointIntervalSequence":
        """0-length DIS at the coordinate space's 5' end."""
        return self._at_index(self.coordinate_end5_index, on_coordinate_strand=True)

    @property
    def coord_end3(self) -> "DisjointIntervalSequence":
        """0-length DIS at the coordinate space's 3' end."""
        return self._at_index(self.coordinate_end3_index, on_coordinate_strand=True)

    # -------------------------------------------------------------------
    # Other properties
    # -------------------------------------------------------------------

    @property
    def coordinate_intervals(self) -> tuple[Interval, ...]:
        """The underlying genomic intervals of the coordinate-space, sorted 5'->3'."""
        return self._coordinate_intervals

    @property
    def coordinate_length(self) -> int:
        """Total length of the coordinate space in bases."""
        return self._coordinate_length

    @property
    def length(self) -> int:
        """Length of the interval in the coordinate space."""
        return abs(self._interval_end3_index - self._interval_end5_index)

    # -------------------------------------------------------------------
    # Private helpers
    # -------------------------------------------------------------------

    def _upstream_index_step(self, on_coordinate_strand: bool | None = None) -> int:
        """Return +1 or -1 indicating the upstream direction in index space.

        Args:
            on_coordinate_strand: Override for which strand to compute the step for.
                Defaults to this interval's on_coordinate_strand.
        """
        if on_coordinate_strand is None:
            on_coordinate_strand = self.on_coordinate_strand
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return -1 if on_coordinate_strand else 1
        # POSITIVE_STRAND_LEFT_TO_RIGHT: effective strand determines direction
        # + strand coord + on_coord_strand → effective +, upstream = -1
        # + strand coord + off_coord_strand → effective -, upstream = +1
        # - strand coord + on_coord_strand → effective -, upstream = +1
        # - strand coord + off_coord_strand → effective +, upstream = -1
        on_positive = (
            (self.coord_transcript_strand == "+" and on_coordinate_strand)
            or (self.coord_transcript_strand == "-" and not on_coordinate_strand)
        )
        return -1 if on_positive else 1

    def _validate_same_coordinate_space(self, other: "DisjointIntervalSequence") -> None:
        """Raise if other does not share the same coordinate space."""
        if not isinstance(other, DisjointIntervalSequence):
            raise TypeError(f"Expected DisjointIntervalSequence, got {type(other).__name__}")
        if self._coordinate_intervals != other._coordinate_intervals:
            raise ValueError("DIS objects must share the same coordinate intervals")

    def _from_indices(self, end5: int, end3: int) -> "DisjointIntervalSequence":
        """Return a new DIS with the same coordinate space but different interval indices."""
        return DisjointIntervalSequence(
            self._coordinate_intervals,
            coord_id=self._coord_metadata.id,
            interval_id=self._interval_metadata.id,
            on_coordinate_strand=self.on_coordinate_strand,
            interval_end5_index=end5,
            interval_end3_index=end3,
        )

    # -------------------------------------------------------------------
    # Interval transform methods
    # -------------------------------------------------------------------

    def shift(self, amount: int) -> "DisjointIntervalSequence":
        """Shift the interval downstream by amount (negative shifts upstream).

        The coordinate space is unchanged. Only the interval indices move.
        """
        downstream_step = -self._upstream_index_step()
        delta = amount * downstream_step
        return self._from_indices(
            self._interval_end5_index + delta,
            self._interval_end3_index + delta,
        )

    def expand(self, upstream: int, dnstream: int | None = None) -> "DisjointIntervalSequence":
        """Expand the interval upstream and/or downstream.

        Negative values contract the interval. Raises ValueError if contraction
        would result in end5 being downstream of end3.

        Args:
            upstream: Bases to expand (or contract if negative) toward the 5' end.
            dnstream: Bases to expand (or contract if negative) toward the 3' end.
                Defaults to upstream (symmetric).
        """
        if dnstream is None:
            dnstream = upstream
        up_step = self._upstream_index_step()
        new_end5 = self._interval_end5_index + upstream * up_step
        new_end3 = self._interval_end3_index - dnstream * up_step
        # Validate end5 is still upstream of or equal to end3
        if (new_end5 - new_end3) * up_step < 0:
            raise ValueError(
                "Invalid expansion: end5 would be downstream of end3 "
                f"(end5={new_end5}, end3={new_end3})"
            )
        return self._from_indices(new_end5, new_end3)

    # -------------------------------------------------------------------
    # Relational methods
    # -------------------------------------------------------------------

    def upstream_of(self, other: "DisjointIntervalSequence") -> bool:
        """True if self is strictly upstream of other (no overlap).

        Requires the same coordinate space and same on_coordinate_strand.
        """
        self._validate_same_coordinate_space(other)
        if self.on_coordinate_strand != other.on_coordinate_strand:
            raise ValueError("Cannot compare: intervals are on different strands")
        if (self.length == 0 and other.length == 0
                and self._interval_end5_index == other._interval_end5_index):
            return False
        up_step = self._upstream_index_step()
        # self's 3' end must be at or upstream of other's 5' end
        return (self._interval_end3_index - other._interval_end5_index) * up_step >= 0

    def dnstream_of(self, other: "DisjointIntervalSequence") -> bool:
        """True if self is strictly downstream of other (no overlap).

        Requires the same coordinate space and same on_coordinate_strand.
        """
        self._validate_same_coordinate_space(other)
        if self.on_coordinate_strand != other.on_coordinate_strand:
            raise ValueError("Cannot compare: intervals are on different strands")
        if (self.length == 0 and other.length == 0
                and self._interval_end5_index == other._interval_end5_index):
            return False
        up_step = self._upstream_index_step()
        # self's 5' end must be at or downstream of other's 3' end
        return (self._interval_end5_index - other._interval_end3_index) * up_step <= 0

    def within(self, other: "DisjointIntervalSequence") -> bool:
        """True if self's interval is contained within other's interval.

        Requires the same coordinate space and same on_coordinate_strand.
        """
        self._validate_same_coordinate_space(other)
        if self.on_coordinate_strand != other.on_coordinate_strand:
            raise ValueError("Cannot compare: intervals are on different strands")
        up_step = self._upstream_index_step()
        # self's end5 must be downstream of or equal to other's end5
        end5_ok = (self._interval_end5_index - other._interval_end5_index) * up_step <= 0
        # self's end3 must be upstream of or equal to other's end3
        end3_ok = (self._interval_end3_index - other._interval_end3_index) * up_step >= 0
        return end5_ok and end3_ok

    # -------------------------------------------------------------------
    # Strand methods
    # -------------------------------------------------------------------

    def is_positive_strand(self) -> bool:
        """True if the interval is on the positive strand."""
        if self.transcript_strand == "+":
            return True
        return False

    def as_positive_strand(self) -> "DisjointIntervalSequence":
        """Return a DIS with the interval on the positive strand. Returns self if already positive."""
        if self.is_positive_strand():
            return self
        return self.as_opposite_strand()

    def as_negative_strand(self) -> "DisjointIntervalSequence":
        """Return a DIS with the interval on the negative strand. Returns self if already negative."""
        if not self.is_positive_strand():
            return self
        return self.as_opposite_strand()

    def as_opposite_strand(self) -> "DisjointIntervalSequence":
        """Return a new DIS with the interval on the opposite strand."""
        return DisjointIntervalSequence(
            self._coordinate_intervals,
            coord_id=self._coord_metadata.id,
            interval_id=self._interval_metadata.id,
            on_coordinate_strand=not self.on_coordinate_strand,
            interval_end5_index=self._interval_end3_index,
            interval_end3_index=self._interval_end5_index,
        )

    # -------------------------------------------------------------------
    # Other methods
    # -------------------------------------------------------------------

    def genomic_span(self) -> Interval:
        """Smallest single genomic Interval representing the DIS, that spans the disjoint genomic coordinate intervals."""
        pass

    def __len__(self) -> int:
        """Return the length of the interval."""
        return self.length

    def __repr__(self) -> str:
        """Return a human-readable representation."""
        return (
            f"DisjointIntervalSequence("
            f"coord_id={self._coord_metadata.id!r}, "
            f"id={self._interval_metadata.id!r}, "
            f"{self.chromosome}:{self.coord_transcript_strand}, "
            f"len={self.length}, "
            f"coord_intervals={len(self._coordinate_intervals)})"
            f"end5={self.interval_end5_index}, "
            f"end3={self.interval_end3_index})"
        )

    def __eq__(self, other: object) -> bool:
        """Equality based on coordinate intervals, metadata, and index values."""
        if not isinstance(other, DisjointIntervalSequence):
            return NotImplemented
        return (
            self._coordinate_intervals == other._coordinate_intervals
            and self._coord_metadata == other._coord_metadata
            and self._interval_metadata == other._interval_metadata
            and self._interval_end5_index == other._interval_end5_index
            and self._interval_end3_index == other._interval_end3_index
        )
