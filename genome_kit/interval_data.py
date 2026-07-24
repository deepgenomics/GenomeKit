from .diseq import DisjointIntervalSequence
from .interval import Interval
from typing import TypeAlias, Sequence

# Interval-like types IntervalData can be backed by and indexed with. Both present
# the same 5'->3' coordinate-space interface (__len__, strand, end5.start, expand,
# as_opposite_strand, contains), so IntervalData stays agnostic to which one it holds.
_INTERVAL_LIKE = (Interval, DisjointIntervalSequence)
IntervalLike: TypeAlias = Interval | DisjointIntervalSequence


class IntervalData:
    r"""Associates a data sequence with a genomic interval for convenient splicing.

    Array indexing via ``[index]`` acts directly on ``data``. Array splicing (via
    ``slice``, :py:class:`~genome_kit.Interval`, or
    :py:class:`~genome_kit.DisjointIntervalSequence`) splices both the ``interval``
    and ``data``; the associated ``data`` is assumed to be unstranded 5'->3', so
    splicing via the opposite strand will reverse the ``data`` (similar to tracks).

    Examples
    --------
    >>> from genome_kit import Interval, IntervalData
    >>> import numpy as np
    >>> interval_len = 500
    >>> position = Interval('chr7', '+', 100000, 100000, 'hg19')
    >>> interval = position.expand(0, interval_len)
    >>> data = np.arange(5 * interval_len * 10).reshape(5, interval_len, 10)
    >>> interval_data = IntervalData(interval, data, axis=1)
    >>>
    >>> # Slicing with __getitem__ is applied to both the interval and data array
    >>> str(interval_data[:, 100:200, 2:6])
    '<<type \'numpy.ndarray\'> indexed by <Interval("chr7", "+", 100100, 100200, "hg19")>>'
    >>>
    >>> # Slicing with interval reproduces the same results
    >>> str(interval_data[position.shift(100).expand(0, 100)])
    '<<type \'numpy.ndarray\'> indexed by <Interval("chr7", "+", 100100, 100200, "hg19")>>'
    >>>
    >>> # Indexing with __getitem__ is applied to the data array
    >>> interval_data[1, 0, 0]
    5000

    """

    def __init__(self, interval: IntervalLike, data, axis=0):
        """Initialize an IntervalData.

        Parameters
        ----------
        interval
            The genomic :py:class:`~genome_kit.Interval` or
            :py:class:`~genome_kit.DisjointIntervalSequence` of interest.
        data
            An ``array_like`` the same length as ``interval`` along ``axis`` and
            ordered 5'->3'. If required, ``axis`` can be used to realign;
            otherwise, :func:`~numpy.rollaxis` can be used to reindex.
        axis
            Axis of the multidimensional ``data`` aligned to ``interval``.
            Defaults to 0.

        Raises
        ------
        IndexError
            If ``axis`` is negative.
        ValueError
            If the length of ``data`` along ``axis`` does not match ``interval``.
        TypeError
            If ``interval`` is neither an :py:class:`~genome_kit.Interval` nor a
            :py:class:`~genome_kit.DisjointIntervalSequence`.
        """

        if not isinstance(interval, _INTERVAL_LIKE):
            raise TypeError(
                "interval must be an Interval or DisjointIntervalSequence, got {}.".format(
                    type(interval).__name__
                )
            )

        if axis < 0:
            raise IndexError("axis {} must be non-negative.".format(axis))
        if (
            axis == 0
            and len(interval) != len(data)
            or axis > 0
            and len(interval) != data.shape[axis]
        ):
            raise ValueError(
                "interval ({}) and data ({}) must be of the same rank.".format(
                    len(interval), len(data)
                )
            )
        self.axis = axis

        self.interval = interval
        self.data = data

    @classmethod
    def from_interval(cls, interval: Interval, data, axis=0) -> "IntervalData":
        """Construct an IntervalData from a single genomic Interval.

        Parameters
        ----------
        interval
            The genomic :py:class:`~genome_kit.Interval` of interest.
        data
            See :py:meth:`__init__`.
        axis
            See :py:meth:`__init__`.
        """
        return cls(interval, data, axis)

    @classmethod
    def from_dis(cls, dis: DisjointIntervalSequence, data, axis=0) -> "IntervalData":
        """Construct an IntervalData from a DisjointIntervalSequence.

        Parameters
        ----------
        dis
            The :py:class:`~genome_kit.DisjointIntervalSequence` whose segment
            defines the genomic interval(s).
        data
            See :py:meth:`__init__`.
        axis
            See :py:meth:`__init__`.
        """
        return cls(dis, data, axis)

    @classmethod
    def from_intervals(
        cls, intervals: Sequence[Interval], data, axis=0
    ) -> "IntervalData":
        """Construct an IntervalData from a sequence of disjoint genomic Intervals.

        The intervals are transformed into a :py:class:`~genome_kit.DisjointIntervalSequence`.

        Parameters
        ----------
        intervals
            Non-overlapping :py:class:`~genome_kit.Interval` objects on the same
            chromosome, strand, and reference genome.
        data
            See :py:meth:`__init__`.
        axis
            See :py:meth:`__init__`.
        """
        return cls(DisjointIntervalSequence.from_intervals(intervals), data, axis)

    def __len__(self) -> int:
        """Return the number of aligned positions, i.e. ``len(interval)``."""
        return len(self.interval)

    @staticmethod
    def _get_interval(interval: IntervalLike, slice_index: slice) -> IntervalLike:
        """Return the sub-interval selected by a slice along the aligned axis.

        A negative step reverses the slice and flips the result to the opposite
        strand, so it stays ordered 5'->3'.

        Parameters
        ----------
        interval
            The interval-like object being sliced.
        slice_index
            The slice applied to the aligned axis.

        Raises
        ------
        KeyError
            If ``slice_index`` has a step other than ``±1``.
        """
        length = len(interval)
        start, stop, step = slice_index.indices(length)
        if abs(step) != 1:
            raise KeyError("discontinuous steps are not supported for intervals.")

        if step > 0:
            upstream = -start
            downstream = stop - length
        else:
            upstream = -(stop + 1)
            downstream = start - length + 1
        interval = interval.expand(upstream, downstream)
        if step < 0:
            interval = interval.as_opposite_strand()
        return interval

    @staticmethod
    def _get_slice(interval: IntervalLike, interval_slice: IntervalLike) -> slice:
        """Return the aligned-axis slice that selects an interval-like key.

        A key on the opposite strand produces a reversed (negative-step) slice.

        Parameters
        ----------
        interval
            The backing interval-like object.
        interval_slice
            The interval-like key to locate within ``interval``.

        Raises
        ------
        IndexError
            If ``interval`` does not contain ``interval_slice``.
        """
        if interval_slice.strand == interval.strand:
            if not interval.contains(interval_slice):
                raise IndexError(
                    "interval {} does not contain {}".format(interval, interval_slice)
                )

            start = abs(interval_slice.end5.start - interval.end5.start)
            stop = start + len(interval_slice)
            step = None
        else:
            opp_strand = interval_slice.as_opposite_strand()
            if not interval.contains(opp_strand):
                raise IndexError(
                    "interval {} does not contain {}".format(interval, interval_slice)
                )

            stop = abs(opp_strand.end5.start - interval.end5.start) - 1
            start = stop + len(opp_strand)
            if stop < 0:
                stop = None
            step = -1
        return slice(start, stop, step)

    def _lift_key(self, key: IntervalLike) -> IntervalLike:
        """Normalize an interval-like key into the backing interval's coordinate space.

        When this IntervalData is backed by a
        :py:class:`~genome_kit.DisjointIntervalSequence` and ``key`` is a genomic
        :py:class:`~genome_kit.Interval`, the key is lifted (via
        :py:meth:`~genome_kit.DisjointIntervalSequence.lift_interval`) so that
        indexing happens in the DIS's flattened coordinate space. All other cases
        — an Interval-backed IntervalData, or a key already in the backing
        coordinate space — pass through unchanged.

        Raises
        ------
        IndexError
            If ``key`` is not fully contained within the backing DIS's segment
            (see :py:meth:`~genome_kit.DisjointIntervalSequence.lift_interval`).
        """
        if isinstance(self.interval, DisjointIntervalSequence) and isinstance(
            key, Interval
        ):
            # lift_interval raises ValueError when `key` is not contained in the
            # segment or is on a mismatched chromosome/reference genome. From the
            # indexing caller's perspective those are all "key not in this data",
            # so surface an IndexError to match the Interval-backed path.
            try:
                lifted = self.interval.lift_interval(key)
            except ValueError as ex:
                raise IndexError(str(ex)) from ex
            if lifted is None:
                raise IndexError(
                    "interval {} does not contain {}".format(self.interval, key)
                )
            return lifted
        return key

    def __getitem__(
        self, item: slice | tuple | int | IntervalLike
    ) -> "IntervalData" | Sequence:
        """Index the data, or splice both the interval and data together.

        An integer (or a tuple whose aligned-axis entry is an integer) indexes
        ``data`` directly and returns the raw array value. A ``slice``, an
        :py:class:`~genome_kit.Interval`, or a
        :py:class:`~genome_kit.DisjointIntervalSequence` slices both ``interval``
        and ``data`` and returns a new :py:class:`IntervalData`. Splicing onto
        the opposite strand reverses ``data``.
        """
        data = self.data
        axis = self.axis

        if isinstance(item, slice):
            if axis > 0:
                raise IndexError(
                    "aligned axis {} requires multidimensional slice.".format(axis)
                )
            return IntervalData(self._get_interval(self.interval, item), data[item])
        elif isinstance(item, tuple):
            index = item[axis]
            if isinstance(index, slice):
                return IntervalData(
                    self._get_interval(self.interval, index), data[item], axis
                )
            return data[item]
        elif isinstance(item, _INTERVAL_LIKE):
            item = self._lift_key(item)
            slice_index = self._get_slice(self.interval, item)
            if axis > 0:
                slices = data.ndim * [slice(None)]
                slices[axis] = slice_index
                slice_index = tuple(slices)
            return IntervalData(item, data[slice_index], axis)
        return data[item]

    def __setitem__(self, key: slice | tuple | IntervalLike, value):
        """Assign into ``data`` in place.

        An :py:class:`~genome_kit.Interval` or
        :py:class:`~genome_kit.DisjointIntervalSequence` key is first resolved to
        the corresponding slice of the aligned axis; a ``slice`` or tuple is
        forwarded to ``data`` directly.
        """
        if isinstance(key, _INTERVAL_LIKE):
            key = self._get_slice(self.interval, self._lift_key(key))
        self.data[key] = value

    def __repr__(self):
        """Return a human-readable representation."""
        return "IntervalData({}, {}, {})".format(
            repr(self.interval), repr(self.data), repr(self.axis)
        )

    def __str__(self):
        """Return the data type together with the interval it is indexed by."""
        return "<{} indexed by <{}>>".format(type(self.data), self.interval)
