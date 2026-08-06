from .diseq import DisjointIntervalSequence
from .interval import Interval
from typing import TypeAlias, Sequence
import numpy as np

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
    '<<class \'numpy.ndarray\'> indexed by <Interval("chr7", "+", 100100, 100200, "hg19")>>'
    >>>
    >>> # Slicing with interval reproduces the same results
    >>> str(interval_data[position.shift(100).expand(0, 100)])
    '<<class \'numpy.ndarray\'> indexed by <Interval("chr7", "+", 100100, 100200, "hg19")>>'
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
                "interval ({}) and data axis ({}) must be of the same length.".format(
                    len(interval), len(data) if axis == 0 else data.shape[axis]
                )
            )
        self._axis = axis

        self._interval = interval
        self._data = data

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

    def __len__(self) -> int:
        """Return the number of aligned positions, i.e. ``len(interval)``."""
        return len(self._interval)

    @property
    def interval(self) -> IntervalLike:
        """Return the backing interval-like object."""
        return self._interval

    @property
    def data(self):
        """Return the backing data array."""
        return self._data

    @property
    def axis(self) -> int:
        """Return the axis of ``data`` aligned to ``interval``."""
        return self._axis

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

    @staticmethod
    def _validate_interval_list(items: list) -> None:
        """Validate a ``list`` key/item of genomic Intervals.

        Every element must be an :py:class:`~genome_kit.Interval`, and no two
        elements may overlap (a shared boundary, e.g. ``[10, 20)`` and
        ``[20, 30)``, is fine).

        Raises
        ------
        TypeError
            If any element is not an :py:class:`~genome_kit.Interval`.
        ValueError
            If any two elements overlap.
        """
        if not all(isinstance(item, Interval) for item in items):
            raise TypeError("list elements must be Interval.")
        if len(items) < 2:
            return
        if items[0].strand == "+":
            sorted_items = sorted(items, key=lambda iv: iv.start)
        else:
            sorted_items = sorted(items, key=lambda iv: -iv.end)
        for a, b in zip(sorted_items, sorted_items[1:]):
            if a.overlaps(b):
                raise ValueError(
                    "list elements must not overlap: {} overlaps {}.".format(a, b)
                )

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
        if isinstance(self._interval, DisjointIntervalSequence) and isinstance(
            key, Interval
        ):
            # lift_interval raises ValueError when `key` is not contained in the
            # segment or is on a mismatched chromosome/reference genome. From the
            # indexing caller's perspective those are all "key not in this data",
            # so surface an IndexError to match the Interval-backed path.
            try:
                lifted = self._interval.lift_interval(key)
            except ValueError as ex:
                raise IndexError(str(ex)) from ex
            if lifted is None:
                raise IndexError(
                    "interval {} does not contain {}".format(self._interval, key)
                )
            return lifted
        return key

    def __getitem__(
        self, item: slice | tuple | int | IntervalLike | list[Interval]
    ) -> "IntervalData" | Sequence:
        """Index the data, or splice both the interval and data together.

        An integer (or a tuple whose aligned-axis entry is an integer) indexes
        ``data`` directly and returns the raw array value. A ``slice``, an
        :py:class:`~genome_kit.Interval`, a
        :py:class:`~genome_kit.DisjointIntervalSequence`, or a ``list`` of
        non-overlapping :py:class:`~genome_kit.Interval` objects slices both
        ``interval`` and ``data`` and returns a new :py:class:`IntervalData`.
        Splicing onto the opposite strand reverses ``data``.
        """
        data = self._data
        axis = self._axis

        if isinstance(item, slice):
            if axis > 0:
                raise IndexError(
                    "aligned axis {} requires multidimensional slice.".format(axis)
                )
            return IntervalData(self._get_interval(self._interval, item), data[item])
        elif isinstance(item, tuple):
            index = item[axis]
            if isinstance(index, slice):
                return IntervalData(
                    self._get_interval(self._interval, index), data[item], axis
                )
            return data[item]
        elif isinstance(item, _INTERVAL_LIKE) or isinstance(item, list):
            return self._get_interval_like(item)
        return data[item]

    def _get_interval_like(self, item: IntervalLike | list[Interval]) -> "IntervalData":
        """Return the sub-``IntervalData`` selected by an interval-like item.

        If ``item`` is a :py:class:`~genome_kit.DisjointIntervalSequence` or ``List``,
        the multi-region data returned from this IntervalData is
        concatenated along the aligned axis, in 5'->3' order.

        Raises
        ------
        TypeError
            If ``item`` is a ``list`` containing an element that is not an
            :py:class:`~genome_kit.Interval`.
        ValueError
            If ``item`` is a ``list`` containing two overlapping Intervals.
        """
        if isinstance(item, list):
            self._validate_interval_list(item)
            dis = DisjointIntervalSequence.from_intervals(item)
            return IntervalData(dis, self._concat_pieces(dis.lower()), self._axis)

        if isinstance(self._interval, Interval) and isinstance(item, DisjointIntervalSequence):
            return IntervalData(item, self._concat_pieces(item.lower()), self._axis)

        lifted = self._lift_key(item)
        slice_index = self._get_slice(self._interval, lifted)
        if self._axis > 0:
            slices = self._data.ndim * [slice(None)]
            slices[self._axis] = slice_index
            slice_index = tuple(slices)
        return IntervalData(lifted, self._data[slice_index], self._axis)

    def _concat_pieces(self, intervals: list[Interval]):
        """Fetch and concatenate the data slices for ``intervals`` (in order)
        along the aligned axis.
        """
        pieces = [self._get_interval_like(iv).data for iv in intervals]
        return pieces[0] if len(pieces) == 1 else np.concatenate(pieces, axis=self._axis)

    # The forms a value assigned via an interval-like key may take. See
    # :py:meth:`_resolve_value_form`.
    _SCALAR = "scalar"
    _BLOCK = "block"
    _POSITION = "position"

    def _data_shape(self) -> tuple:
        return self._data.shape if hasattr(self._data, "shape") else np.shape(self._data)

    @staticmethod
    def _value_shape(value) -> tuple | None:
        """Return ``value``'s shape, or ``None`` if it has no well-defined one
        (e.g. a ragged nested ``list``).
        """
        if hasattr(value, "shape"):
            return value.shape
        try:
            return np.shape(value)
        except ValueError:
            return None

    def _resolve_value_form(self, value, total: int) -> str:
        """Classify ``value`` against the ``total`` aligned positions selected by
        an interval-like key.

        Exactly three forms are accepted, chosen so that which axis of ``value``
        is the aligned one is never ambiguous — ``value``'s rank pins it down:

        ``"block"``
            ``value`` has ``data``'s shape with the aligned axis resized to
            ``total``, i.e. the shape :py:meth:`__getitem__` returns for the same
            key. Supplies one value per selected position.
        ``"position"``
            ``value`` has ``data``'s shape with the aligned axis removed. Describes
            a single position, and is broadcast to every selected one.
        ``"scalar"``
            ``value`` has no shape at all, and is broadcast to everything.

        A value that would otherwise reach the aligned axis only through NumPy's
        right-aligned broadcasting is rejected, since the axis it lands on depends
        on its rank rather than on the caller's intent.

        Raises
        ------
        ValueError
            If ``value`` matches none of the three forms.
        """
        shape = self._value_shape(value)
        if shape == ():
            return self._SCALAR
        data_shape = self._data_shape()
        block = data_shape[: self._axis] + (total,) + data_shape[self._axis + 1 :]
        position = data_shape[: self._axis] + data_shape[self._axis + 1 :]
        if shape == block:
            return self._BLOCK
        if shape == position:
            return self._POSITION
        raise ValueError(
            f"value shape {shape if shape is not None else 'unknown'} does not fit "
            f"{total} aligned position(s): expected {block} (one value per position), "
            f"{position} (one value broadcast to every position), or a scalar."
        )

    def _normalize_to_backing_coord(self, key: IntervalLike | list[Interval]) -> list[IntervalLike]:
        """Resolve an interval-like key into the region(s) of the aligned axis it
        selects, each normalized into the backing coordinate space and ordered
        5'->3'.

        A ``list`` key is normalized through a
        :py:class:`~genome_kit.DisjointIntervalSequence`, so its elements come back
        merged where adjacent and in 5'->3' order

        Raises
        ------
        TypeError
            If ``key`` is a ``list`` containing an element that is not an
            :py:class:`~genome_kit.Interval`.
        ValueError
            If ``key`` is an empty ``list``, or a ``list`` containing two
            overlapping Intervals.
        IndexError
            If ``key`` is not contained within the backing interval.
        """
        if isinstance(key, list):
            self._validate_interval_list(key)
            key = DisjointIntervalSequence.from_intervals(key)
            return [self._lift_key(iv) for iv in key.lower()]
        if isinstance(key, DisjointIntervalSequence) and isinstance(self._interval, Interval):
            return [self._lift_key(iv) for iv in key.lower()]
        return [self._lift_key(key)]

    def _set_interval_like(self, key: IntervalLike | list[Interval], value) -> None:
        """Assign ``value`` into ``data`` at the aligned-axis region(s) selected by
        an interval-like key.

        A ``"block"`` value is split along the aligned axis, one chunk per selected
        region in 5'->3' order; the other forms are broadcast to every region.

        Raises
        ------
        TypeError
            If ``key`` is a ``list`` containing an element that is not an
            :py:class:`~genome_kit.Interval`.
        ValueError
            If ``key`` is an empty ``list`` or contains overlapping Intervals, or if
            ``value`` matches none of the accepted forms.
        IndexError
            If ``key`` is not contained within the backing interval.
        """
        keys = self._normalize_to_backing_coord(key)
        lengths = [len(interval) for interval in keys]
        form = self._resolve_value_form(value, sum(lengths))

        if form == self._BLOCK:
            offset = 0
            for interval, length in zip(keys, lengths):
                self._set_slice(interval, self._slice_axis(value, offset, offset + length))
                offset += length
            return
        if form == self._POSITION:
            # Re-insert the aligned axis so NumPy broadcasts `value` along it;
            # left as-is it would right-align onto the trailing axes instead.
            value = np.expand_dims(value, self._axis)
        for interval in keys:
            self._set_slice(interval, value)

    def _set_slice(self, key: IntervalLike, value) -> None:
        """Assign ``value`` into ``data`` at the aligned-axis slice selected by a
        single ``key`` already normalized into the backing coordinate space.
        """
        slice_index = self._get_slice(self._interval, key)
        if self._axis > 0:
            slices = len(self._data_shape()) * [slice(None)]
            slices[self._axis] = slice_index
            slice_index = tuple(slices)
        self._data[slice_index] = value

    def _slice_axis(self, value, start: int, stop: int):
        """Return ``value`` sliced to ``[start:stop)`` along the aligned axis
        of this IntervalData.
        """
        if hasattr(value, "shape"):
            slices = value.ndim * [slice(None)]
            slices[self._axis] = slice(start, stop)
            return value[tuple(slices)]
        return value[start:stop]

    def __setitem__(self, key: slice | tuple | IntervalLike | list[Interval], value):
        """Assign into ``data`` in place.

        An :py:class:`~genome_kit.Interval`, a
        :py:class:`~genome_kit.DisjointIntervalSequence`, or a ``list`` of
        non-overlapping :py:class:`~genome_kit.Interval` objects is resolved
        to the corresponding region(s) of the aligned axis;
        a ``slice`` or tuple is forwarded to ``data`` directly.

        For an interval-like ``key``, ``value`` must take one of three forms
        (see :py:meth:`_resolve_value_form`): ``data``'s shape with the aligned
        axis resized to the key's length, supplying one value per selected
        position; ``data``'s shape with the aligned axis removed, describing a
        single position to broadcast to every selected one; or a scalar. Shapes
        that would reach the aligned axis only via NumPy's right-aligned
        broadcasting — including ones padded with length-1 axes — are rejected.

        Raises
        ------
        IndexError
            If ``key`` is a plain ``slice`` and the aligned axis is not axis 0,
            or if an interval-like ``key`` is not contained within ``interval``.
        TypeError
            If ``key`` is a ``list`` containing an element that is not an
            :py:class:`~genome_kit.Interval`.
        ValueError
            If ``key`` is an empty ``list`` or a ``list`` containing two
            overlapping Intervals, or if ``value`` matches none of the accepted
            forms.
        """
        if isinstance(key, _INTERVAL_LIKE) or isinstance(key, list):
            self._set_interval_like(key, value)
            return
        if isinstance(key, slice) and self._axis > 0:
            raise IndexError(
                "aligned axis {} requires multidimensional slice.".format(self._axis)
            )
        self._data[key] = value

    def __repr__(self):
        """Return a human-readable representation."""
        return "IntervalData({}, {}, {})".format(
            repr(self._interval), repr(self._data), repr(self._axis)
        )

    def __str__(self):
        """Return the data type together with the interval it is indexed by."""
        return "<{} indexed by <{}>>".format(type(self._data), self._interval)
