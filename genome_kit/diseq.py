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
class _IntervalMetadata:
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
    - An **interval** within that coordinate space, defined by a 5' and 3' index.
      The interval may lie on the same, or opposite, strand as the coordinate space.

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
        interval_name: str | None = None,
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
        interval_name : :py:class:`str` or None
            Optional name for the interval.
        on_coordinate_strand : :py:class:`bool`
            Whether the interval is on the same strand as the coordinate
            intervals.
        start : :py:class:`int` or None
            start index of the interval in the coordinate space. Defaults to 0
        end : :py:class:`int` or None
            end index of the interval in the coordinate space. Defaults to the length
            of the coordinate space.

        Raises
        ------
        ValueError
            If intervals are empty, inconsistent, overlapping, or if start
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
        self._interval_metadata = _IntervalMetadata(
            name=interval_name,
            on_coordinate_strand=on_coordinate_strand,
        )

        # Default interval start/end to span the full coordinate
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
        interval_name: str | None = None,
    ) -> "DisjointIntervalSequence":
        """Construct a DIS from a sequence of Intervals
        (or :py:class:`~genome_kit.Exon`/:py:class:`~genome_kit.Cds`/:py:class:`~genome_kit.Utr` objects).

        If elements have an ``.interval`` attribute (e.g. Exon, Cds, Utr),
        the plain Interval is extracted automatically.

        Parameters
        ----------
        intervals : Sequence[:py:class:`~genome_kit.Interval`]
            Sequence of Interval or annotation objects with ``.interval``.
        coord_name : :py:class:`str` or None
            Optional name for the coordinate space.
        interval_name : :py:class:`str` or None
            Optional name for the interval.

        Returns
        -------
        :py:class:`DisjointIntervalSequence`
        """
        # Extract .interval if items are Exon/Cds/Utr
        coord_intervals = []
        for iv in intervals:
            if type(iv) is not Interval and hasattr(iv, "interval"):
                coord_intervals.append(iv.interval)
            else:
                coord_intervals.append(iv)
        return cls(coord_intervals, coord_name=coord_name, interval_name=interval_name)

    @classmethod
    def from_transcript(
        cls,
        transcript: Transcript,
        *,
        region: Literal["exons", "cds", "utr5", "utr3"] = "exons",
        coord_name: str | None = None,
        interval_name: str | None = None,
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
        interval_name : :py:class:`str` or None
            Optional name for the interval. Defaults to ``transcript.id``.

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
        if interval_name is None:
            interval_name = transcript.id

        return cls(coord_intervals, coord_name=coord_name, interval_name=interval_name)

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
    def coord_transcript_strand(self) -> Literal["+", "-"]:
        """Strand of the coordinate intervals (the transcript strand)."""
        return self._coord_metadata.transcript_strand

    @property
    def name(self) -> str | None:
        """Name of the interval, or None."""
        return self._interval_metadata.name

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

    @property
    def coordinate_end5_index(self) -> int:
        """5' index of the coordinate space."""
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return 0
        if self.coord_transcript_strand == "+":
            return 0
        return self.coordinate_length

    @property
    def coordinate_end3_index(self) -> int:
        """3' index of the coordinate space."""
        if self._index_direction == IndexDirection.TRANSCRIPT_FIVE_TO_THREE:
            return self.coordinate_length
        if self.coord_transcript_strand == "+":
            return self.coordinate_length
        return 0

    @property
    def end5_index(self) -> int:
        """5' index of the interval."""
        if self._upstream_index_step() == -1:
            return self._start
        return self._end

    @property
    def end3_index(self) -> int:
        """3' index of the interval."""
        if self._upstream_index_step() == -1:
            return self._end
        return self._start

    @property
    def start(self) -> int:
        """Start index of the interval in the coordinate space. Not necessarily the 5' end."""
        return self._start

    @property
    def end(self) -> int:
        """End index of the interval in the coordinate space. Not necessarily the 3' end."""
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
        """0-length DIS at the interval's 5' end."""
        return self._at_index(
            self.end5_index, on_coordinate_strand=self.on_coordinate_strand
        )

    @property
    def end3(self) -> "DisjointIntervalSequence":
        """0-length DIS at the interval's 3' end."""
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
        """Length of the interval on the coordinate space."""
        return self.end - self.start

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
        return -1 if self.transcript_strand == "+" else 1

    def __len__(self) -> int:
        """Return the length of the interval."""
        return self.length

    def __repr__(self) -> str:
        """Return a human-readable representation."""
        return (
            f"DisjointIntervalSequence("
            f"coord_name={self._coord_metadata.name!r}, "
            f"name={self._interval_metadata.name!r}, "
            f"{self.chromosome}:{self.coord_transcript_strand}, "
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
        try:
            # If refg mismatch, ValueError is raised
            self._coordinate_intervals == other._coordinate_intervals
        except ValueError:
            return False
        return (
            self._coord_metadata == other._coord_metadata
            and self._interval_metadata == other._interval_metadata
            and self._start == other._start
            and self._end == other._end
            and self._coordinate_intervals == other._coordinate_intervals
        )
