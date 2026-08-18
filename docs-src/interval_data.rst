.. _interval_data:

-------------------------------
Interval Data
-------------------------------

Overview
========

An :py:class:`~genome_kit.IntervalData` is an object that associates a NumPy array with
a genomic interval. Conceptually, it is similar to a
:py:class:`~genome_kit.GenomeTrack`, but in-memory (doesn't require a file) and
scoped to a single interval-like object.

More specifically, an :py:class:`~genome_kit.IntervalData` binds three things:

- ``interval`` — the backing coordinate object, either an
  :py:class:`~genome_kit.Interval` or a
  :py:class:`~genome_kit.DisjointIntervalSequence`.
- ``data`` — an array whose aligned axis has the same length as ``interval``.
- ``axis`` — which axis of ``data`` is aligned to ``interval`` (default ``0``).

The data is assumed to be **unstranded and ordered 5'→3'** along the interval,
much like a :ref:`track <tracks>`. Indexing then falls into two categories:

- **Array indexing** (an integer, or a tuple that indexes the aligned axis with
  an integer) acts directly on ``data`` and returns the raw array element.
- **Splicing** (a ``slice``, an :py:class:`~genome_kit.Interval`, or a
  :py:class:`~genome_kit.DisjointIntervalSequence`) slices *both* the interval
  and the data, returning a new :py:class:`~genome_kit.IntervalData` whose
  ``interval`` and ``data`` still line up. Splicing onto the opposite strand
  reverses the data, so index 0 always corresponds to the 5' end of the result.

.. code-block:: python

    >>> import genome_kit as gk
    >>> import numpy as np
    >>> interval = gk.Interval("chr7", "+", 100000, 100500, "hg19")
    >>> data = np.arange(len(interval))
    >>> interval_data = gk.IntervalData(interval, data)
    >>> len(interval_data)
    500

Construction
============

An :py:class:`~genome_kit.IntervalData` can be backed two ways: by an
:py:class:`~genome_kit.Interval` for a contiguous genomic region, or by a
:py:class:`~genome_kit.DisjointIntervalSequence` for a discontiguous one (e.g. a
transcript's exons). Either way, the aligned axis of ``data`` must have the same
length as the backing object — for a DIS that is the length of the segment
(``len(dis)``)

.. code-block:: python

    >>> # backed by an Interval: one value per base of a contiguous region
    >>> interval = gk.Interval("chr7", "+", 100000, 100500, "hg19")
    >>> data = np.arange(len(interval))          # length 500
    >>> interval_data = gk.IntervalData(interval, data)
    >>>
    >>> # backed by a DIS: one value per base of the spliced sequence
    >>> from genome_kit.diseq import DisjointIntervalSequence
    >>> dis = DisjointIntervalSequence.from_transcript(transcript)
    >>> data = np.arange(len(dis))
    >>> interval_data = gk.IntervalData.from_dis(dis, data)

The constructor takes either kind of backing object;
:py:meth:`~genome_kit.IntervalData.from_interval` and
:py:meth:`~genome_kit.IntervalData.from_dis` are explicit aliases that behave
identically to it.

The aligned axis
~~~~~~~~~~~~~~~~~

Data is rarely one-dimensional. The ``axis`` argument tells
:py:class:`~genome_kit.IntervalData` which axis is aligned to the interval

.. code-block:: python

    >>> interval = gk.Interval("chr7", "+", 100000, 100500, "hg19")
    >>> data = np.zeros((5, len(interval), 10))   # length is on axis 1
    >>> interval_data = gk.IntervalData(interval, data, axis=1)
    >>> len(interval_data)
    500

Array Indexing
==============

Indexing that resolves to a concrete position on the aligned axis returns the
underlying array data directly

.. code-block:: python

    >>> data = np.arange(30).reshape(3, 10)
    >>> interval = gk.Interval("chr1", "+", 100, 103, "hg19")
    >>> interval_data = gk.IntervalData(interval, data)
    >>> interval_data[0]           # first position along the interval
    array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

Iterating an :py:class:`~genome_kit.IntervalData` yields these per-position
values in order, so ``[x for x in interval_data]`` walks the data 5'→3'.

A tuple whose entry on the aligned axis is an integer likewise returns raw array
data — the whole tuple is forwarded to the underlying array

.. code-block:: python

    >>> data = np.arange(30).reshape(5, 3, 2)
    >>> interval = gk.Interval("chr1", "+", 100, 103, "hg19")
    >>> interval_data = gk.IntervalData(interval, data, axis=1)
    >>> interval_data[1, 0, 0]     # aligned axis (1) indexed by an int -> raw value
    6

Splicing
========

Splicing produces a new :py:class:`~genome_kit.IntervalData` in which the
``interval`` and the ``data`` have been sliced together, so they remain aligned.
There are three ways to splice.

By slice
~~~~~~~~

A Python ``slice`` selects a sub-range along the aligned axis. The backing
interval is narrowed to the matching genomic sub-range

.. code-block:: python

    >>> interval = gk.Interval("chr1", "+", 100, 103, "hg19")
    >>> data = np.arange(30).reshape(3, 10)
    >>> interval_data = gk.IntervalData(interval, data)
    >>> sub = interval_data[1:]            # drop the first position
    >>> sub.interval
    Interval("chr1", "+", 101, 103, "hg19")
    >>> len(sub)
    2

A negative step reverses the slice **and** flips the strand of the interval, so
the data and the interval stay consistent — position 0 of the result is still
its 5' end

.. code-block:: python

    >>> rev = interval_data[::-1]
    >>> rev.interval.strand
    '-'

Only steps of ``±1`` are supported. A discontinuous step (e.g. ``[::2]``) raises
``KeyError``, since no contiguous interval can describe every other base.

When ``axis > 0``, a bare slice is ambiguous (it does not say what to do with the
other axes), so ``interval_data[1:]`` raises ``IndexError``. Use a tuple that
slices the aligned axis explicitly instead

.. code-block:: python

    >>> data = np.arange(30).reshape(5, 3, 2)
    >>> interval_data = gk.IntervalData(interval, data, axis=1)
    >>> sub = interval_data[:, 1:, ...]    # slice the aligned axis (1)
    >>> sub.interval
    Interval("chr1", "+", 101, 103, "hg19")
    >>> np.array_equal(interval_data[:, 1:, ...].data, data[:, 1:, ...])
    True  # A slice is forwarded to the underlying array

By Interval or DIS
~~~~~~~~~~~~~~~~~~~

The most convenient form is to index with a genomic key. The key selects the
portion of the data corresponding to that genomic region, and the returned
:py:class:`~genome_kit.IntervalData` is indexed by the key itself

.. code-block:: python

    >>> base = gk.Interval("chr1", "+", 100, 103, "hg19")
    >>> data = np.arange(30).reshape(3, 10)
    >>> interval_data = gk.IntervalData(base, data)
    >>> key = base.expand(-1, 0)           # drop one base at the 5' end
    >>> sub = interval_data[key]
    >>> sub.interval == key
    True
    >>> np.array_equal([x for x in sub], data[1:])
    True

The key must be contained within the backing interval; otherwise ``IndexError``
is raised.

Indexing with the interval on the **opposite strand** reverses the data, so the
result is read 5'→3' along the requested strand

.. code-block:: python

    >>> rev = interval_data[base.as_opposite_strand()]
    >>> rev.interval.strand
    '-'
    >>> np.array_equal([x for x in rev], data[::-1])
    True

A zero-length key (e.g. ``base.end5``) yields an empty result, and a
full-length key round-trips to the original data.

Lifting genomic keys onto a DIS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the :py:class:`~genome_kit.IntervalData` is backed by a
:py:class:`~genome_kit.DisjointIntervalSequence`, you can index it with a plain
genomic :py:class:`~genome_kit.Interval`. The key is first *lifted* into the DIS
coordinate space via
:py:meth:`~genome_kit.DisjointIntervalSequence.lift_interval`, so a genomic region
is translated into spliced-coordinate offsets before the data is sliced

.. code-block:: python

    >>> dis = DisjointIntervalSequence.from_transcript(transcript)
    >>> data = np.arange(len(dis))
    >>> interval_data = gk.IntervalData.from_dis(dis, data)
    >>> exon = transcript.exons[0].interval
    >>> sub = interval_data[exon]          # genomic key lifted into DIS space

The genomic key must be fully contained within the DIS segment, otherwise
an ``IndexError`` is raised.

It is also possible to index a DIS-backed :py:class:`~genome_kit.IntervalData` with a
:py:class:`~genome_kit.DisjointIntervalSequence` key, but they must share the same
coordinate intervals.

Assigning Data
==============

Assignment mirrors indexing. A key that is a ``slice`` or a tuple is forwarded
directly to the underlying array

.. code-block:: python

    >>> interval_data[0] = -interval_data[0]

An :py:class:`~genome_kit.Interval` or
:py:class:`~genome_kit.DisjointIntervalSequence` key is first resolved to the
corresponding slice of the aligned axis (lifting the key onto the backing DIS if
necessary), then assigned into ``data``

.. code-block:: python

    >>> base = gk.Interval("chr1", "+", 100, 103, "hg19")
    >>> data = np.arange(30).reshape(3, 10)
    >>> interval_data = gk.IntervalData(base, data.copy())
    >>> interval_data[base.expand(-1, 0)] = 0    # zero out all but the first position
    >>> np.array_equal(interval_data.data[0], data[0])
    True
    >>> bool(np.all(interval_data.data[1:] == 0))
    True

A ``list`` of non-overlapping :py:class:`~genome_kit.Interval` objects assigns
across several regions at once. The list is normalized 5'->3' first, so the order
it is given in does not matter.

When assigning via an Interval-like or ``List`` object, ``value`` must
take one of three shapes:

- ``data``'s shape with the aligned axis resized to the key's length — one value
  per selected position, i.e. the shape of ``interval_data[key].data``
- ``data``'s shape with the aligned axis removed — a single position's worth of
  data, broadcast to every position the key selects.
- a scalar, broadcast to everything.

.. code-block:: python

    >>> exons = [gk.Interval("chr1", "+", 100, 102, "hg19"),
    ...          gk.Interval("chr1", "+", 104, 106, "hg19")]
    >>> base = gk.Interval("chr1", "+", 100, 106, "hg19")
    >>> interval_data = gk.IntervalData(base, np.zeros((6, 10), dtype=int))
    >>> interval_data[exons] = np.arange(40).reshape(4, 10)   # one value per position
    >>> interval_data[exons] = np.arange(10)                  # one row, to all four
    >>> interval_data[exons] = 0                              # scalar

Any other shape raises ``ValueError`` to prevent silent bugs due to ambiguous assignment.

Assignment mutates the existing ``data`` array in place; it does not create a new
:py:class:`~genome_kit.IntervalData`.

Attributes and Representation
=============================

An :py:class:`~genome_kit.IntervalData` exposes its components directly:

- ``interval`` — the backing :py:class:`~genome_kit.Interval` or
  :py:class:`~genome_kit.DisjointIntervalSequence`.
- ``data`` — the underlying array.
- ``axis`` — the aligned axis.

See Also
========

- :ref:`diseq` — the coordinate system used when an
  :py:class:`~genome_kit.IntervalData` is backed by a
  :py:class:`~genome_kit.DisjointIntervalSequence`.
- :ref:`tracks <tracks>` — similar to an :py:class:`~genome_kit.IntervalData`,
  but on-disk and for a whole genome.
