from typing import NamedTuple
import numpy as np


class KeyLayout(NamedTuple):
    """Where one tracked axis of an array goes when the array is indexed by a key.

    Attributes
    ----------
    expanded_key
        The key with its ``Ellipsis`` expanded and trailing ``:`` added, so every
        axis of the array has an entry (plus any ``None``).
    axis_key_index
        The index into ``expanded_key`` of the entry for the tracked axis.
    new_axis
        The axis of the result the tracked axis becomes, or None if the tracked
        axis is indexed by an integer (and so dropped).
    """

    expanded_key: list
    axis_key_index: int
    new_axis: int | None


# What ``key_layout`` accepts, for its error message.
BASIC_INDICES = "integers, slices, None and Ellipsis"


def _is_advanced(entry) -> bool:
    # Types used for numpy advanced indexing
    return isinstance(entry, (bool, np.bool_, list, tuple, np.ndarray))


def _is_integer(entry) -> bool:
    # bool is a subclass of int, but NumPy treats a bool index as a mask.
    return isinstance(entry, (int, np.integer)) and not isinstance(entry, (bool, np.bool_))


def key_layout(key: tuple, ndim: int, axis: int, supported: str = BASIC_INDICES) -> KeyLayout:
    """Work out where ``axis`` of an ``ndim``-dimensional array goes when the array
    is indexed by ``key``, following NumPy's basic indexing rules.

    Only basic indexing (integers, slices, ``None`` and ``Ellipsis``) is supported.
    ``supported`` describes the accepted keys in the ``TypeError`` message, for
    callers that accept more than this function does.

    Raises
    ------
    IndexError
        If ``key`` has more than one ``Ellipsis`` or indexes too many axes.
    NotImplementedError
        If an entry of ``key`` is an advanced index (an index array, a nested
        sequence or a mask).
    TypeError
        If an entry of ``key`` is any other unsupported type.
    """
    for entry in key:
        if _is_advanced(entry):
            raise NotImplementedError(
                "advanced indexing (index arrays, sequences and masks) is not "
                "supported, got {!r}.".format(entry)
            )
        if not (entry is None or entry is Ellipsis or isinstance(entry, slice) or _is_integer(entry)):
            raise TypeError(
                "unsupported index {!r}: only {} are supported.".format(entry, supported)
            )
    if sum(entry is Ellipsis for entry in key) > 1:
        raise IndexError("an index can only have a single ellipsis ('...')")
    consumed = sum(entry is not None and entry is not Ellipsis for entry in key)
    if consumed > ndim:
        raise IndexError(
            "too many indices for array: array is {}-dimensional, but {} were indexed".format(
                ndim, consumed
            )
        )

    # Pad with ``:`` (at the Ellipsis, if any) so every axis has an entry.
    fill = (ndim - consumed) * [slice(None)]
    if Ellipsis in key:
        at = key.index(Ellipsis)
        expanded_key = list(key[:at]) + fill + list(key[at + 1 :])
    else:
        expanded_key = list(key) + fill

    # Find the entry for `axis`, counting the result dimensions before it:
    # slices and Nones keep a dimension, integers drop one.
    dim = 0
    new_axis = 0
    for position, entry in enumerate(expanded_key):
        if entry is None:
            new_axis += 1
            continue
        if dim == axis:
            break
        if isinstance(entry, slice):
            new_axis += 1
        dim += 1
    axis_key_index = position

    if not isinstance(expanded_key[axis_key_index], slice):
        return KeyLayout(expanded_key, axis_key_index, None)
    return KeyLayout(expanded_key, axis_key_index, new_axis)
