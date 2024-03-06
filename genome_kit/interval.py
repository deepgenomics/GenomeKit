# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

from operator import attrgetter
from warnings import warn

from . import _cxx
from ._cxx_util import mock, mock_result


@_cxx.register
class Interval(_cxx.Interval):
    """A genomic interval.

    An Interval represents a contiguous stranded genomic interval, identified
    by a chromosome, strand, start position, end position, and reference genome.
    """

    __slots__ = ()  # <--- DO NOT EXTEND BASE CLASS SLOTS

    # __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # __del__(self):  pass   # <--- DO NOT IMPLEMENT
    # __getattribute__       # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # __setattribute__       # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION

    def __init__(self, chromosome, strand, start, end, reference_genome, anchor=None, anchor_offset=0):
        """Initialize an interval from a DNA0 position.

        Parameters
        ----------
        chromosome : :py:class:`str`
            The chromosome name ('chr1', ...).

        strand : :py:class:`chr`
            The strand ("+" or "-").

        start : :py:class:`int`
            The start position of the interval (inclusive, 0-based).

        end : :py:class:`int`
            The end position of the interval (exclusive, 0-based).

        reference_genome : :py:class:`str` | :py:class:`~genome_kit.Genome`
            The reference genome in which the interval is defined,
            for example "hg19". If a genome object is given, its reference
            genome will be used.

        anchor : :py:class:`int` | :py:class:`str` | :py:class:`~genome_kit.Interval` | None, optional
            The anchor is only relevant when extracting DNA from a
            `VariantGenome` and specifies how length-changing variants are
            applied. Conceptually, the anchor defines a fixed point when the
            sequence changes length. See (TODO: link) for the details.
            The anchor can be any position in the interval. If anchor is
            `None`, then length-changing variants will change the length of the
            interval. Anchor can be a string: `5p` or `3p` sets the anchor to
            the respective end of the interval, respecting the sense strand.
            `start` or `end` set the anchor to the start or end of the interval
            on the forward strand. When anchor is `center` then the anchor is
            set to the middle position of the interval (if the interval has odd
            length, the position is rounded towards the 5p end of the interval,
            respecting the sense strand). If anchor is an integer or an empty
            interval, the anchor is set to that position.

        anchor_offset : :py:class:`int`, optional
            The anchor offset is a shift that's applied to the anchor only on
            the variant genome. The purpose is to align a position on the
            reference genome with a position *inside* an insertion. The anchor
            itself would only allow to align with the beginning or end of the
            insertion.
        """

        if isinstance(anchor, str):
            if anchor == 'start':
                anchor = start
            elif anchor == 'end':
                anchor = end
            elif anchor == '5p':
                anchor = start if strand == '+' else end
            elif anchor == '3p':
                anchor = end if strand == '+' else start
            elif anchor == 'center':
                anchor = start + (end - start) // 2 if strand == '+' else \
                         end   - (end - start) // 2  # noqa
            else:
                raise ValueError("Unrecognized string '{}' for anchor".format(anchor))

        # Deliberately skip overhead of super(Interval, self) wrapper
        _cxx.Interval.__init__(self, chromosome, strand, start, end, reference_genome, anchor, anchor_offset)

    @mock
    @property
    def chromosome(self):  # pragma: no cover
        """The chromosome name.

        A shorthand property `chrom` is also available.

        Returns
        -------
        :py:class:`str`
            The name of the chromosome this interval is defined on  (e.g. "chr1").
        """
        return mock_result(str)

    @mock
    @property
    def strand(self):  # pragma: no cover
        """The strand ("+" or "-").

        Returns
        -------
        :py:class:`chr`
            The strand this interval is defined on.
        """
        return mock_result(chr)

    @mock
    @property
    def start(self):  # pragma: no cover
        """The integer position of the start (inclusive, 0-based).

        Returns
        -------
        :py:class:`int`
            The start position of this interval within the reference genome.
        """
        return mock_result(int)

    @mock
    @property
    def end(self):  # pragma: no cover
        """The integer position of the end (exclusive, 0-based).

        Returns
        -------
        :py:class:`int`
            The end position of the this interval within the reference genome.
        """
        return mock_result(int)

    @mock
    @property
    def end5(self):  # pragma: no cover
        """The 5p end of this interval, which is a length-0 interval.

        Not to be confused with the `end` attribute, which


        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The length-0 interval at the 5p end of the interval, i.e. the gap preceding
            the most upstream base in the interval.
        """
        return mock_result(Interval)

    @mock
    @property
    def end3(self):  # pragma: no cover
        """The 3p end of this interval, which is a length-0 interval.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The length-0 interval at the 3p end of the interval, i.e. the gap succeeding
            the most downstream base in the interval.
        """
        return mock_result(Interval)

    @mock
    @property
    def reference_genome(self):  # pragma: no cover
        """The reference genome (assembly name in UCSC format).

        A shorthand property `refg` is also available.

        Returns
        -------
        :py:class:`str`
            The name of the reference genome this interval is defined in.
        """
        return mock_result(str)

    @mock
    @property
    def anchor(self):  # pragma: no cover
        """The anchor position.

        Returns
        -------
        :py:class:`int` | :py:data:`None`
            The anchor position associated with this interval, if any.
        """
        return mock_result(int)

    @mock
    @property
    def anchor_offset(self):  # pragma: no cover
        """The anchor offset.

        Returns
        -------
        :py:class:`int`
            The anchor offset associated with this interval, if any.
        """
        return mock_result(int)

    @mock
    def shift(self, amount):  # pragma: no cover
        """The interval shifted upstream/downstream.

        Parameters
        ----------
        amount : :py:class:`int`
            The amount by which to shift the start and end positions.
            A positive amount shifts downstream according to the current strand,
            and a negative amount shifts the position upstream.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The shifted interval.
        """
        return mock_result(Interval)

    @mock
    def expand(self, upstream, dnstream=None):  # pragma: no cover
        """The interval expanded upstream/downstream.

        Can also be called with a single argument, in which case
        upstream and dnstream take on the same value.

        Parameters
        ----------
        upstream : :py:class:`int`
            Number of upstream positions to add to the interval.

        dnstream : :py:class:`int`
            Number of downstream positions to add to the interval.
            If not specified, defaults to `upstream` (i.e. symmetric).

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            An interval containing `upstream` + `downstream` + 1 positions.
        """
        return mock_result(Interval)

    @mock
    def intersect(self, other):  # pragma: no cover
        """The interval intersecting both intervals.

        Parameters
        ----------
        other : :class:`~genome_kit.Interval`
            interval to intersect with this interval.

        Returns
        -------
        :class:`~genome_kit.Interval` | :data:`None`
            The intersecting interval (possibly empty) or `None` if the
            intervals are disjoint
        """
        return mock_result(Interval)

    def subtract(self, other: Interval) -> Sequence[Interval]:
        """A list of intervals, representing this interval with the its
        intersection of another interval removed.

        Parameters
        ----------
        other
            interval to subtract its intersection from this interval.

        Returns
        -------
        Sequence[Interval]
            The intervals after subtraction. If the intervals were disjoint,
            returns [self].
        """
        intersection = self.intersect(other)
        if not intersection:
            return [self]
        result = []
        if self.start < intersection.start:
            result.append(Interval(self.chrom, self.strand, self.start, intersection.start, self.refg))
        if self.end > intersection.end:
            result.append(Interval(self.chrom, self.strand, intersection.end, self.end, self.refg))
        return result

    @mock
    def is_positive_strand(self):  # pragma: no cover
        """Test if this interval is on the positive strand.

        Returns
        -------
        :py:class:`bool`
            True if `strand` is equal to '+'
        """
        return mock_result(bool)

    @mock
    def as_positive_strand(self):  # pragma: no cover
        """The same interval on the positive strand.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The same interval, but on the positive ('+') strand.
        """
        return mock_result(Interval)

    @mock
    def as_negative_strand(self):  # pragma: no cover
        """The same interval on the negative strand.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The same interval, but on the negative ('-') strand.
        """
        return mock_result(Interval)

    @mock
    def as_opposite_strand(self):  # pragma: no cover
        """The same interval on the opposite strand.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The same interval, but on the opposite strand.
        """
        return mock_result(Interval)

    @mock
    def upstream_of(self, other):  # pragma: no cover
        """Test if this interval is strictly upstream of `other`.

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval`
            Any coordinate or interval subclass defined on the same
            reference genome and strand as this interval.

        Returns
        -------
        :py:class:`bool`
            True if this interval is strictly upstream of `other`, with no overlap.
        """
        return mock_result(bool)

    @mock
    def dnstream_of(self, other):  # pragma: no cover
        """Test if this interval is strictly downstream of `other`.

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval`
            Any coordinate or interval subclass defined on the same
            reference genome and strand as this interval.

        Returns
        -------
        :py:class:`bool`
            True if this interval is strictly downstream of `other`, with no overlap.
        """
        return mock_result(bool)

    @mock
    def contains(self, other):  # pragma: no cover
        """Test if `other` lies within the extents of this interval.

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval`
            Any coordinate or interval subclass defined on the same
            reference genome and strand as this interval.

        Returns
        -------
        :py:class:`bool`
            True if `other` lies entirely within the extents of this interval.
        """
        return mock_result(bool)

    @mock
    def within(self, other):  # pragma: no cover
        """Test if this interval lies within the extents of `other`.

        Use of ``interval in other`` expressions is also available.

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval`
            Any coordinate or interval subclass defined on the same
            reference genome and strand as this interval.

        Returns
        -------
        :py:class:`bool`
            True if this interval lies within the extents of `other`.
        """
        return mock_result(bool)

    @mock
    def overlaps(self, other):  # pragma: no cover
        """Test if this interval overlaps the extents of `other`.

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval`
            Any coordinate or interval subclass defined on the same
            reference genome and strand as this interval.

        Returns
        -------
        :py:class:`bool`
            True if this interval overlaps the extents of `other`.
        """
        return mock_result(bool)

    @staticmethod
    def from_dna0_coord(chromosome, strand, position, reference_genome):
        """Create an interval spanning a single DNA0 position.

        Equivalent to
        ``Interval(chromosome, strand, position, position+1, reference_genome)``.

        Parameters
        ----------
        chromosome : :py:class:`str`
            The chromosome name (e.g. "chr1").

        strand : :py:class:`int`
            The strand ("+" or "-").

        position : :py:class:`int`
            The position the interval should span (0-based).

        reference_genome : :py:class:`str` | :py:class:`~genome_kit.Genome`
            The reference genome in which the interval is defined,
            for example "hg19". If a genome object is given, its reference
            genome will be used.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The resulting length-1 interval spanning `position`.
        """
        return Interval(chromosome, strand, position, position + 1, reference_genome)

    @staticmethod
    def from_dna0(chromosome, strand, start, end, reference_genome):
        """Create an interval from DNA0 positions.

        Equivalent to calling :py:class:`~genome_kit.Interval`.

        Parameters
        ----------
        chromosome : :py:class:`str`
            The chromosome name (e.g. "chr1").

        strand : :py:class:`int`
            The strand ("+" or "-").

        start : :py:class:`int`
            The start position of the interval (inclusive, 0-based).

        end : :py:class:`int`
            The end position of the interval (exclusive, 0-based).

        reference_genome : :py:class:`str` | :py:class:`~genome_kit.Genome`
            The reference genome in which the interval is defined,
            for example "hg19". If a genome object is given, its reference
            genome will be used.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The resulting interval.
        """
        return Interval(chromosome, strand, start, end, reference_genome)

    @staticmethod
    def from_dna1(chromosome, strand, start, end, reference_genome):
        """Create an interval from DNA1 positions.

        Positions are automatically converted to DNA0 as (start-1, end).

        Parameters
        ----------
        chromosome : :py:class:`str`
            The chromosome name (e.g. "chr1").

        strand : :py:class:`int`
            The strand ("+" or "-").

        start : :py:class:`int`
            The start position of the interval (inclusive, 1-based).

        end : :py:class:`int`
            The end position of the interval (inclusive, 1-based).

        reference_genome : :py:class:`str` | :py:class:`~genome_kit.Genome`
            The reference genome in which the interval is defined,
            for example "hg19". If a genome object is given, its reference
            genome will be used.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The resulting interval.
        """

        return Interval(chromosome, strand, start - 1, end, reference_genome)

    @staticmethod
    def from_rna1(chromosome, start, end, reference_genome):
        """Create an interval from RNA1 positions.

        The start/end positions are automatically converted to DNA0 as
        `min(abs(start),abs(end))-1` and `max(abs(start), abs(end)))` respectively.

        Parameters
        ----------
        chromosome : :py:class:`str`
            The chromosome name (e.g. "chr1").

        strand : :py:class:`int`
            The strand ("+" or "-").

        start : :py:class:`int`
            The start position of the interval (signed, inclusive, 1-based).

        end : :py:class:`int`
            The end position of the interval (signed, inclusive, 1-based).

        reference_genome : :py:class:`str` | :py:class:`~genome_kit.Genome`
            The reference genome in which the interval is defined,
            for example "hg19". If a genome object is given, its reference
            genome will be used.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The resulting interval.
        """

        strand = '+' if start > 0 else '-'

        # verify input
        if end:
            if not start <= end:
                raise ValueError("Start position needs to be <= end (%d, %d)" % (start, end))

            if (start >= 0) != (end >= 0):
                raise ValueError("Start and end must have the same sign: " "(%d, %d)" % (start, end))

        # check if it's single position or interval
        start, end = min(abs(start), abs(end)) - 1, max(abs(start), abs(end))

        return Interval(chromosome, strand, start, end, reference_genome)

    @staticmethod
    def spanning(interval1, interval2):
        """Create the smallest interval spanning the given pair of existing intervals,
        or raise an exception if given intervals not on the same chromosome, strand,
        and reference genome. The given intervals are not required to be overlapping,
        and may both be of zero length.

        Example::

            >>> a = Interval("chr1", "+", 100, 150, "hg19")
            >>> b = Interval("chr1", "+", 200, 250, "hg19")
            >>> Interval.spanning(a, b)
            Interval("chr1", "+", 100, 250, "hg19")

        Notes
        -----
            This method is not anchor aware, since that requires a list of variants for context.
            See :func:`~genome_kit.Variant.spanning` for that use case.

        Parameters
        ----------
        interval1 : :py:class:`~genome_kit.Interval`
            an interval

        interval2 : :py:class:`~genome_kit.Interval`
            another interval


        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The resulting spanning interval.
        """

        interval1.check_same_chromosome(interval2)
        if interval1.strand != interval2.strand:
            warn(f"Spanning intervals on opposite strands, got: {interval1}, {interval2}")
        if any(x.anchor is not None for x in [interval1, interval2]):
            raise ValueError(
                "Spanning anchored intervals (should use Variant.spanning), got:"
                f" {interval1}, {interval2}"
            )

        start = min(interval1.start, interval2.start)
        end = max(interval1.end, interval2.end)

        return Interval(interval1.chrom, interval1.strand, start, end, interval1.reference_genome)

    @property
    def midpoint(self):
        """Identify the midpoint of an interval.

        Example::

            >>> Interval("chr1", "+", 100, 150, "hg19").midpoint
            Interval("chr1", "+", 125, 125, "hg19")
            >>> Interval("chr1", "+", 100, 151, "hg19").midpoint
            Interval("chr1", "+", 125, 126, "hg19")

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanning the midpoint of the interval.
        """
        return self.expand(-(len(self) // 2))

    def distance(self, other, method='midpoint'):
        """The distance this interval and another.

        Example::

            >>> interval_a = Interval("chr1", "+", 100, 150, "hg19")
            >>> interval_b = Interval("chr1", "+", 200, 300, "hg19")
            >>> interval_c = Interval("chr1", "+", 400, 500, "hg19")
            >>> interval_a.distance(interval_b)
            125.0
            >>> interval_a.distance([interval_b, interval_c])
            [125.0, 325.0]

        Parameters
        ----------
        other : :py:class:`~genome_kit.Interval` or :py:class:`~list`
            an interval or list of intervals
        method : :py:class:`str`
            determines which (point) attribute of this and the other interval to
            use when computing the distance.
            Must be one of `midpoint`, `end5`, or `end3`.

        Returns
        -------
        :py:class:`float` or :py:class:`~list`
            Distance to the interval. If other is of type Interval, return distance as a float;
            otherwise, a list of distances.
        """
        def midpoint_as_float(interval):
            x = interval.midpoint
            return x.start + len(x) / 2

        if method == 'midpoint':
            getter = midpoint_as_float
        elif method in ['end5', 'end3']:
            getter = attrgetter(method + '.start')
        else:
            raise ValueError("`method` must be in ['midpoint', 'end5', 'end3']")

        def validated(interval):
            self.check_same_chromosome(interval)
            if self.strand != interval.strand:
                warn(
                    f"Distancing intervals on opposite strands, got: {self}, {interval}"
                )
            if any(x.anchor is not None for x in [self, interval]):
                raise ValueError(
                    f"Distancing anchored intervals is not supported, got: {self},"
                    f" {interval}"
                )
            return interval

        if isinstance(other, Interval):
            return float(abs(getter(self) - getter(validated(other))))
        start = getter(self)
        return [float(abs(start - getter(validated(x)))) for x in other]

    def as_dna0(self):
        """Returns (start, end) positions using DNA0 convention.

        Returns
        -------
        :py:class:`tuple`
            The positions ``(start, end)`` using DNA0 convention.
        """

        return self.start, self.end

    def as_dna1(self):
        """Returns (start, end) positions using DNA1 convention.

        Returns
        -------
        :py:class:`tuple`
            The positions ``(start, end)`` using DNA1 convention.
        """

        if self.start == self.end:
            raise ValueError("Empty intervals cannot be represented using DNA1 convention.")

        return self.start + 1, self.end

    def as_rna1(self):
        """Returns (start, end) positions using RNA1 convention.

        Returns
        -------
        :py:class:`tuple`
            The positions ``(start, end)`` using RNA1 convention.
        """

        if self.start == self.end:
            raise ValueError("Empty intervals cannot be represented using RNA1 convention.")

        if self.is_positive_strand():  # Forward strand
            return self.start + 1, self.end
        else:  # Reverse strand
            return -self.end, -(self.start + 1)

    def as_ucsc(self):
        """Returns this interval as a UCSC browser coordinate.

        Returns
        -------
        :py:class:`str`
            A string in UCSC browser convention, e.g. ``"chr5:500201-5002010"``.
        """

        return "{}:{}-{}".format(self.chromosome, self.start + 1, self.end)

    def with_anchor(self, anchor, anchor_offset=0):
        """Returns an anchored version of this interval.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            An anchored version of this interval.
            The previous anchor, if any, is ignored.
        """

        return Interval(self.chromosome, self.strand, self.start, self.end, self.reference_genome, anchor,
                        anchor_offset)

    def check_same_chromosome(self, other: Interval) -> None:
        if self.chromosome != other.chromosome:
            raise ValueError(
                f"Both intervals must be on the same chromosome, got: {self}, {other}"
            )
        if self.reference_genome != other.reference_genome:
            raise ValueError(
                f"Both intervals must be on the same reference genome, got: {self},"
                f" {other}"
            )

    def __repr__(self):
        anchor_args = ""
        if self.anchor is not None:
            anchor_args += ", " + str(self.anchor)
            if self.anchor_offset:
                anchor_args += ", " + str(self.anchor_offset)
        return 'Interval("{}", "{}", {}, {}, "{}"{})'.format(self.chrom, self.strand, self.start, self.end, self.refg,
                                                             anchor_args)

    def __str__(self):
        return self.__repr__()
