from bisect import bisect_left
from ._indexing import key_layout
from numpy.lib.mixins import NDArrayOperatorsMixin
import numpy as np


class DisjointView(NDArrayOperatorsMixin):
    """A writable view over several disjoint regions of one array.

    NumPy cannot express a view over gapped regions, so this wraps the per-region
    views and presents them as a read-only concatenation along ``axis``. The
    concatenation is built on-demand when requested, so it is kept in sync with the
    ``parent_data``. Assignments are written back through the regions into
    ``parent_data``. Reads (and any other operation) are forwarded to the
    read-only concatenation.

    Basic indexing (ints, slices, ``None`` and ``Ellipsis``) returns another
    ``DisjointView`` over the selected parts of the regions, so chained writes
    such as ``view[1:4][0] = 1`` reach ``parent_data``. An int on ``axis``
    selects from a single region and returns a plain NumPy view of it. Advanced
    indexing (index arrays and masks) raises ``NotImplementedError``, for reads
    and writes.

    Attempting to write to a returned read-only concatenation raises a ValueError.

    Parameters
    ----------
    parent_data
        The original array the regions are views of.
    data
        The per-region views of ``parent_data``, in the order they are
        concatenated.
    axis
        The axis the regions are concatenated along.
    """

    def __init__(self, parent_data, data: list[np.ndarray], axis: int):
        self.parent_data = parent_data
        self.data = data
        self.axis = axis

    def _concatenate(self) -> np.ndarray:
        """Return a fresh, writable concatenation of the regions."""
        return np.concatenate(self.data, axis=self.axis)

    @property
    def concatenated_data(self) -> np.ndarray:
        """Return a read-only concatenation of the regions."""
        concatenated = self._concatenate()
        concatenated.flags.writeable = False
        return concatenated

    def _write_back(self, concatenated: np.ndarray) -> None:
        """Write ``concatenated`` into the regions, and so into ``parent_data``."""
        offset = 0
        for piece in self.data:
            length = piece.shape[self.axis]
            slices = concatenated.ndim * [slice(None)]
            slices[self.axis] = slice(offset, offset + length)
            piece[...] = concatenated[tuple(slices)]
            offset += length

    def __setitem__(self, key, value):
        # Only to reject unsupported keys, as NumPy would accept them.
        key_layout(key if isinstance(key, tuple) else (key,), self.data[0].ndim, self.axis)
        concatenated = self._concatenate()
        concatenated[key] = value
        self._write_back(concatenated)

    def __getitem__(self, key):
        key = key if isinstance(key, tuple) else (key,)
        layout = key_layout(key, self.data[0].ndim, self.axis)
        expanded_key = layout.expanded_key
        axis_key_index = layout.axis_key_index

        axis_length = sum(piece.shape[self.axis] for piece in self.data)
        axis_key = expanded_key[axis_key_index]
        if not isinstance(axis_key, slice):
            index = int(axis_key)
            if not -axis_length <= index < axis_length:
                raise IndexError(
                    "index {} is out of bounds for axis {} with size {}".format(index, self.axis, axis_length)
                )
            index %= axis_length
            for piece in self.data:
                if index < piece.shape[self.axis]:
                    expanded_key[axis_key_index] = index
                    return piece[tuple(expanded_key)]
                index -= piece.shape[self.axis]

        pieces = []
        for piece, local_slice in self._split(range(*axis_key.indices(axis_length))):
            expanded_key[axis_key_index] = local_slice
            pieces.append(piece[tuple(expanded_key)])
        return DisjointView(self.parent_data, pieces, layout.new_axis)

    def _split(self, selected: range) -> list:
        """Split positions ``selected`` along ``axis`` into ``(piece, local slice)``
        pairs, in the order ``selected`` visits them."""
        ascending_range = selected if selected.step > 0 else selected[::-1]
        regions = []
        offset = 0
        for piece in self.data:
            axis_length = piece.shape[self.axis]
            piece_start = bisect_left(ascending_range, offset)
            piece_end = bisect_left(ascending_range, offset + axis_length)
            hit = ascending_range[piece_start:piece_end]
            if hit:
                regions.append((piece, range(hit.start - offset, hit.stop - offset, hit.step)))
            offset += axis_length

        if selected.step < 0:
            regions = [(piece, local_slice[::-1]) for piece, local_slice in reversed(regions)]
        if not regions:
            # No regions selected on the aligned axis (e.g. [0:0]). Preserve numpy
            # behaviour which would preserve the shape of the other axes,
            # but with zero length on the aligned axis.
            regions = [(self.data[0], range(0))]
        return [
            (piece, slice(local_slice.start, local_slice.stop if local_slice.stop >= 0 else None, local_slice.step))
            for piece, local_slice in regions
        ]

    def __bool__(self) -> bool:
        return bool(self.concatenated_data)

    def __getattr__(self, name):
        # Private and protocol names (e.g. ``__deepcopy__``, ``__array_interface__``)
        # are not forwarded: copy/pickle probe them before ``__init__`` has run,
        # and NumPy would read a temporary concatenation's buffer directly.
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self.concatenated_data, name)

    def __len__(self) -> int:
        return len(self.concatenated_data)

    def __iter__(self):
        return iter(self.concatenated_data)

    def __array__(self, dtype=None, copy=None):
        if copy is False:
            raise ValueError("a DisjointView cannot be converted to an array without a copy.")
        concatenated = self._concatenate()
        if dtype is not None:
            concatenated = concatenated.astype(dtype, copy=False)
        if not copy:
            concatenated.flags.writeable = False
        return concatenated

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kwargs):
        target = inputs[0] if method == "at" else None
        inputs = tuple(x._concatenate() if isinstance(x, DisjointView) else x for x in inputs)
        if isinstance(target, DisjointView):
            # ``ufunc.at`` modifies its first operand in place instead of via ``out``.
            getattr(ufunc, method)(*inputs, **kwargs)
            target._write_back(inputs[0])
            return None
        if out is None:
            return getattr(ufunc, method)(*inputs, **kwargs)

        # Compute into writable concatenations, then write them back.
        outputs = tuple(x._concatenate() if isinstance(x, DisjointView) else x for x in out)
        getattr(ufunc, method)(*inputs, out=outputs, **kwargs)
        for original, output in zip(out, outputs):
            if isinstance(original, DisjointView):
                original._write_back(output)
        return out[0] if len(out) == 1 else out

    def __repr__(self):
        return "DisjointView({}, axis={})".format(repr(self.concatenated_data), self.axis)
