.. _diseq:

-------------------------------
Disjoint Interval Sequences
-------------------------------

Motivation
==========

When working with transcripts, it is often necessary to operate on the exonic
sequence with introns (or other intervening regions) removed. The remaining exons
form a disjoint set of genomic :py:class:`~genome_kit.Interval` objects.
Operations such as indexing into the spliced sequence or defining sub-ranges
across exon boundaries require a coordinate system that accounts for the gaps
between these intervals.

The :py:class:`~genome_kit.diseq.DisjointIntervalSequence` (DIS) class provides
this coordinate system. While it builds on :py:class:`~genome_kit.Interval`, a DIS
is conceptually distinct in several ways:

- **Spliced coordinate space.** Positions in a DIS are offsets into the
  concatenated exonic (or other) sequence, not genomic coordinates.
- **5'→3' index direction.** DIS indices always increase from 5' to 3' with
  respect to the transcript. Index 0 corresponds to the transcript's 5' end
  regardless of genomic strand. This contrasts with
  :py:class:`~genome_kit.Interval`, where ``start < end`` always holds in
  genomic coordinates, so on the ``-`` strand ``start`` is the 3' end.
- **Same-strand / opposite-strand semantics.** Because a DIS models spliced
  RNA rather than raw DNA, the concept of ``+``/``-`` strand is replaced by
  ``on_coordinate_strand`` (same strand as the transcript) versus opposite
  strand. The underlying genomic strand is accessible via ``coord_strand``,
  but intervals within the DIS are described relative to the coordinate
  space rather than in absolute genomic terms.

Overview
========

A :py:class:`~genome_kit.diseq.DisjointIntervalSequence` (DIS) represents
a flattened coordinate system over a sequence of disjoint genomic intervals.
For example, the exons of a transcript form a disjoint set of genomic intervals that,
when concatenated, represent the spliced RNA sequence.

A DIS has two aspects:

- A **coordinate space**: the underlying genomic
  :py:class:`~genome_kit.Interval` objects (e.g. exons) that define the
  flattened index system. These intervals are sorted 5'→3' and must not overlap.
- An **interval**: a sub-range within that coordinate space, defined by
  a start and end index, where start <= end.

The following examples illustrate how the coordinate space and interval interact,
using both diagrams and code.

Consider a transcript on the + strand with the following genomic layout:
::
    Genomic Coordinates:  153 154 155 156 157 158 159 160 161 162 163 164 165 166 167
    DNA Sequence:       |  A   T   G   C   C   G   C   A   T   G   C   C   G   C  |
                        | |<------->| |<------->| |<--->| |<----------->| |<--->| |
                        5'   Exon1      Intron1    Exon2      Intron2      Exon3  3'

Extracting only the exons yields the following disjoint intervals:
::
    Genomic Coordinates:  153 154 155 159 160 165 166 167
    DNA Sequence:       |  A   T   G   C   A   G   C  |
                        | |<------->| |<--->| |<--->| |
                        5'   Exon1     Exon2   Exon3  3'

These exon intervals can be represented as :py:class:`~genome_kit.Interval` objects::

    >>> from genome_kit import Interval
    >>> from genome_kit.diseq import DisjointIntervalSequence
    >>> exon1 = Interval("chr1", "+", 153, 156, "hg38")
    >>> exon2 = Interval("chr1", "+", 159, 161, "hg38")
    >>> exon3 = Interval("chr1", "+", 165, 167, "hg38")

To define an interval spanning the full exonic sequence (from the start of Exon1 to
the end of Exon3), the exon intervals are first converted into a DIS coordinate space
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:       |  A   T   G   C   A   G   C  |
                        | |<------->| |<--->| |<--->| |
                        5'   Exon1     Exon2   Exon3  3'

The default interval spans the entire coordinate space
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          A   T   G   C   A   G   C
                           |<--------------------->|
                          end5      Interval      end3
    Start Index:     0
    End Index:       7

The interval spans the full length of the coordinate space, with a start index of 0
and an end index of 7::

    >>> dis = DisjointIntervalSequence.from_intervals(
    ...     [exon1, exon2, exon3], coord_name="tx_example"
    ... )
    >>> dis.start
    0
    >>> dis.end
    7
    >>> dis.on_coordinate_strand
    True

.. note::
    The disjoint interval follows the convention of
    :py:class:`~genome_kit.Interval` where intervals are half-open
    (the end index is exclusive).

A DIS can also represent an interval on the strand opposite the coordinate space.
This is useful for modeling the complementary sequence or a binding partner.

Starting from the coordinate space defined above
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:       |  A   T   G   C   A   G   C  |
                        |  |<----->|   |<->|   |<->|  |
                        5'   Exon1     Exon2   Exon3  3'

The opposite strand shares the same DIS coordinate indices
::
                        5'      Positive strand          3'
    DIS Coordinates:    |  0   1   2   3   4   5   6   7 |
    DNA Sequence (+):   |  A   T   G   C   A   G   C     |
    -----------------------------------------------------
    DNA Sequence (-):   |  T   A   C   G   T   C   G     |
    DIS Coordinates:    |  0   1   2   3   4   5   6   7 |
                        3'      Negative Strand          5'

The DIS coordinate indices are identical on both strands. To obtain the complement
of a given interval, the same start and end indices apply; only the
``on_coordinate_strand`` flag changes. The following shows the full-length interval
on the opposite strand
::
                        5'      Coordinate Strand        3'
    DIS Coordinates:    |  0   1   2   3   4   5   6   7 |
    DNA Sequence (+):   |  A   T   G   C   A   G   C     |
    -----------------------------------------------------
    DNA Sequence (-):      T   A   C   G   T   C   G
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                  Opposite Strand
                           |<--------------------->|
                          end3      Interval      end5
    Start Index:     0
    End Index:       7
    On Coordinate Strand: False

The ``on_coordinate_strand`` flag distinguishes same-strand from opposite-strand
intervals, since the start and end indices alone do not encode strand information::

    >>> dis_opp = DisjointIntervalSequence(
    ...     [exon1, exon2, exon3],
    ...     coord_name="tx_example",
    ...     on_coordinate_strand=False,
    ... )
    >>> dis_opp.on_coordinate_strand
    False
    >>> dis_opp.end5_index
    7
    >>> dis_opp.end3_index
    0

.. note::
    The preceding examples used + strand coordinate intervals. When the coordinate intervals
    lie on the **negative** strand, the DIS behaves differently from :py:class:`~genome_kit.Interval`:
    in one key aspect. Indices still increase in the 5'→3' direction of the transcript.::

        >>> iv = Interval("chr1", "-", 0, 100, genome.refg)
        >>> assert iv.end5.start > iv.end3.start
        >>> dis = DisjointIntervalSequence.from_intervals([iv])
        >>> assert dis.end5.start < dis.end3.start

Consider a transcript on the negative strand:
::
                        3'  Exon3      Intron2    Exon2      Intron1      Exon1   5'
                        | |<------->| |<------->| |<--->| |<----------->| |<--->| |
    DNA Sequence (-):   |  G   T   C   A   G   T   C   A   G   T   C   A   G   T  |
    Genomic Coordinates:  153 154 155 156 157 158 159 160 161 162 163 164 165 166 167
                                            Negative Strand
Extracting only the exons:
::
                        3'   Exon3     Exon2   Exon1  5'
                        | |<------->| |<--->| |<--->| |
    DNA Sequence (-):   |  G   T   C   C   A   G   T  |
    Genomic Coordinates:  153 154 155 159 160 165 166 167
                                    Negative Strand

Converting these exons into a DIS coordinate space:
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:       |  T   G   A   C   C   T   G  |
                        |  |<->|   |<->|   |<----->|  |
                        5' Exon1   Exon2     Exon3    3'

The sequence appears reversed relative to genomic coordinates because the DIS
coordinate space is oriented 5'→3' with respect to the transcript, regardless of
genomic strand. Index 0 always corresponds to the transcript's 5' end, and the
largest index to the 3' end::

    >>> neg_exon1 = Interval("chr1", "-", 165, 167, "hg38")
    >>> neg_exon2 = Interval("chr1", "-", 159, 161, "hg38")
    >>> neg_exon3 = Interval("chr1", "-", 153, 156, "hg38")
    >>> dis_neg = DisjointIntervalSequence.from_intervals(
    ...     [neg_exon1, neg_exon2, neg_exon3],
    ...     coord_name="tx_neg_example",
    ... )
    >>> dis_neg.coord_strand
    '-'
    >>> dis_neg.coordinate_length
    7

A full-length interval on the coordinate strand
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
                           |<--------------------->|
                          end5      Interval      end3
    Start Index:     0
    End Index:       7
    On Coordinate Strand: True

Despite creating the DIS from the negative strand, the full-length interval on the
coordinate strand is identical to the + strand example. When working with DIS
objects, strand is expressed only as "same strand" or "opposite strand"::

    >>> dis_neg.start
    0
    >>> dis_neg.end
    7
    >>> dis_neg.on_coordinate_strand
    True

The same coordinate space with an opposite-strand interval
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence (-):      T   G   A   C   C   T   G
    -----------------------------------------------------
    DNA Sequence (+):      A   C   T   G   G   A   C
    DIS Coordinates:       0   1   2   3   4   5   6   7
                           |<--------------------->|
                          end3      Interval      end5
    Start Index:     0
    End Index:       7
    On Coordinate Strand: False

::

    >>> dis_neg_opp = DisjointIntervalSequence(
    ...     [neg_exon1, neg_exon2, neg_exon3],
    ...     coord_name="tx_neg_example",
    ...     on_coordinate_strand=False,
    ... )
    >>> dis_neg_opp.on_coordinate_strand
    False
    >>> dis_neg_opp.end5_index
    7
    >>> dis_neg_opp.end3_index
    0

Construction
============

From a Transcript
~~~~~~~~~~~~~~~~~

The most common way to create a DIS is from a
:py:class:`~genome_kit.Transcript`::

    >>> from genome_kit import Genome
    >>> from genome_kit.diseq import DisjointIntervalSequence
    >>> genome = Genome("gencode.v29")
    >>> transcript = genome.transcripts[100]
    >>> dis = DisjointIntervalSequence.from_transcript(transcript)

By default, the coordinate space is built from the transcript's exons.
It's also possible to specify a region to use CDS or UTR intervals::

    >>> dis_cds = DisjointIntervalSequence.from_transcript(transcript, region="cds")
    >>> dis_utr5 = DisjointIntervalSequence.from_transcript(transcript, region="utr5")
    >>> dis_utr3 = DisjointIntervalSequence.from_transcript(transcript, region="utr3")

The ``coord_name`` and ``interval_name`` default to ``transcript.id`` but can
be overridden::

    >>> dis = DisjointIntervalSequence.from_transcript(
    ...     transcript, coord_name="my_coord", interval_name="my_interval")

From Intervals
~~~~~~~~~~~~~~

A DIS can be constructed from any sequence of
:py:class:`~genome_kit.Interval` objects::

    >>> from genome_kit import Interval
    >>> exon_intervals = [e.interval for e in transcript.exons]
    >>> dis = DisjointIntervalSequence.from_intervals(exon_intervals, coord_name="my_coord")

The intervals must all share the same chromosome, strand, reference
genome, and must not overlap. They are automatically sorted 5'→3'.

Coordinate Space
================

The coordinate space is defined by the underlying genomic intervals, which
are accessible as a tuple::

    >>> dis.coordinate_intervals
    (Interval("chr1", "+", 100, 200, "hg38"), Interval("chr1", "+", 300, 450, "hg38"))
    >>> dis.coordinate_length
    250

Metadata about the coordinate space is available through properties::

    >>> dis.chromosome
    'chr1'
    >>> dis.coord_strand
    '+'
    >>> dis.reference_genome
    'hg38'
    >>> dis.coord_name
    'ENST00000...'

Interval Start and End
======================

The interval within the coordinate space is defined by ``start`` and ``end``
indices, following the same half-open convention as :py:class:`~genome_kit.Interval`
(``start <= end`` always)::

    >>> dis.start
    0
    >>> dis.end
    250
    >>> dis.length
    250
    >>> len(dis)
    250

By default, the interval spans the full coordinate space (``start=0``,
``end=coordinate_length``). Indices can extend beyond ``[0, coordinate_length]``, but
the DNA sequence returned by ``genome.dna()`` will be N-padded.

End5 and End3
~~~~~~~~~~~~~

The ``end5_index`` and ``end3_index`` properties give the 5' and 3' positions
of the interval. These are derived from ``start`` and ``end`` based on the
interval's strand.

When ``on_coordinate_strand`` is ``True``, ``end5_index`` equals ``start`` and
``end3_index`` equals ``end``::

    Start Index:     1
    End Index:       6
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
                               |<------------->|
    -----------------------------------------------------
    DNA Sequence:          A   C   T   G   G   A   C
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                  Opposite Strand

::

    >>> dis = DisjointIntervalSequence.from_transcript(transcript)
    >>> dis.end5_index   # same as start when on coordinate strand
    1
    >>> dis.end3_index   # same as end when on coordinate strand
    6

When ``on_coordinate_strand`` is ``False``, the mapping reverses:
``end5_index`` equals ``end`` and ``end3_index`` equals ``start``::

    Start Index:     3
    End Index:       7
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
    -----------------------------------------------------
    DNA Sequence:          A   C   T   G   G   A   C
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                       |<--------->|
                                  Opposite Strand

::

    >>> opp = dis.as_opposite_strand()
    >>> opp.end5_index   # same as end when off coordinate strand
    7
    >>> opp.end3_index   # same as start when off coordinate strand
    3

Boundary Properties
~~~~~~~~~~~~~~~~~~~

Zero-length DIS objects at the interval and coordinate boundaries are
available as properties::

    >>> dis.end5        # 0-length DIS at the interval's 5' boundary
    >>> dis.end3        # 0-length DIS at the interval's 3' boundary
    >>> dis.coord_end5  # 0-length DIS at the coordinate space's 5' boundary
    >>> dis.coord_end3  # 0-length DIS at the coordinate space's 3' boundary


Strand Methods
==============

A DIS interval can sit on either 'virtual' strand independently of the coordinate
intervals. The ``on_coordinate_strand`` property indicates whether the
interval is on the same strand as the coordinate intervals::
    On Coordinate Strand: True
    Start Index:     1
    End Index:       6
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence (+):      A   T   C   C   G   A   C
                               |<------------->|
    -----------------------------------------------------
    DNA Sequence (-):      T   A   G   G   C   T   G
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                  Opposite Strand

    >>> dis.on_coordinate_strand
    True
    >>> dis.is_same_strand()
    True
    >>> dis.is_positive_strand()
    True

``is_same_strand()`` tests whether the interval is on the coordinate
strand. ``is_positive_strand()`` tests the effective genomic strand
(accounting for both ``coord_strand`` and ``on_coordinate_strand``).

Five methods change the interval's strand. All preserve ``start``,
``end``, and the coordinate intervals. These methods return a DIS on the requested
strand, instead of modifying the existing DIS in-place.

``as_opposite_strand()`` sets ``on_coordinate_strand`` to ``False``::

    Before as_opposite_strand() (on_coordinate_strand=True):
    Start Index:     1
    End Index:       6
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence (+):      T   A   A   C   C   C   T
                               |<------------->|
    -----------------------------------------------------
    DNA Sequence (-):      A   T   T   G   G   G   A
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                  Opposite Strand

    After as_opposite_strand() (on_coordinate_strand=False):
    Start Index:     1
    End Index:       6
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence (+):      T   A   A   C   C   C   T
    -----------------------------------------------------
    DNA Sequence (-):      A   T   T   G   G   G   A
    DIS Coordinates:       0   1   2   3   4   5   6   7
                               |<------------->|
                                  Opposite Strand

::

    >>> dis.on_coordinate_strand
    True
    >>> opposite = dis.as_opposite_strand()
    >>> opposite.on_coordinate_strand
    False
    >>> opposite.start == dis.start   # start/end unchanged
    True

``as_same_strand()`` sets ``on_coordinate_strand`` to ``True``::

    >>> dis.on_coordinate_strand
    False
    >>> dis.is_same_strand()
    False
    >>> same_strand_dis = dis.as_same_strand()
    >>> same_strand_dis.is_same_strand()
    True

``flip_strand()`` toggles ``on_coordinate_strand``::

    >>> dis.on_coordinate_strand
    True
    >>> flipped = dis.flip_strand()
    >>> flipped.on_coordinate_strand
    False
    >>> flipped.flip_strand().on_coordinate_strand
    True

The ``as_positive_strand()`` and ``as_negative_strand()`` methods return a DIS with
the interval on the effective genomic strand::

    >>> dis.coord_strand
    '+'
    >>> dis.on_coordinate_strand
    True
    >>> dis.strand
    '+'
    >>> neg_dis = dis.as_negative_strand()
    >>> neg_dis.strand
    '-'
    >>> pos_dis = neg_dis.as_positive_strand()
    >>> pos_dis.strand
    '+'
    >>> pos_dis.coord_strand == dis.coor_strand == '+'
    True

.. note::

    Strand methods only affect the interval layer. The coordinate
    intervals always remain unchanged.

Shifting and Expanding
======================

Both ``shift`` and ``expand`` return a **new** DIS with modified interval
indices. The coordinate space is always unchanged.

shift
~~~~~

``shift(amount)`` moves the interval downstream by ``amount`` bases.
A negative value shifts upstream. The interval length is preserved.

On the coordinate strand, downstream means increasing indices::

    Before shift(1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                               |<--------->|
                              end5        end3

    After shift(1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                   |<--------->|
                                  end5        end3

On the opposite strand, "downstream" is the reverse direction in index
space, so ``shift(1)`` moves the interval toward *lower* indices::

    Before shift(1) (on_coordinate_strand=False):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                               |<--------->|
                              end3        end5

    After shift(1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                           |<--------->|
                          end3        end5

::

    >>> dis.start, dis.end
    (30, 150)
    >>> shifted = dis.shift(10)
    >>> shifted.start, shifted.end
    (40, 160)
    >>> shifted.coordinate_intervals == dis.coordinate_intervals
    True

    >>> # Negative values shift upstream
    >>> dis.shift(-10).start, dis.shift(-10).end
    (20, 140)

    >>> # On the opposite strand, downstream reverses in index space
    >>> opp = dis.as_opposite_strand()
    >>> opp.start, opp.end
    (30, 150)
    >>> shifted_opp = opp.shift(10)
    >>> shifted_opp.start, shifted_opp.end
    (20, 140)

.. note::

    ``shift`` can move the interval beyond the coordinate space bounds
    (``start < 0`` or ``end > coordinate_length``).

expand
~~~~~~

``expand(upstream, dnstream)`` grows (or shrinks) the interval toward
its 5' and 3' ends. When ``dnstream`` is omitted the expansion is
symmetric::

    Before expand(1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                               |<--------->|
                              end5        end3

    After expand(1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                           |<----------------->|
                          end5                end3

Negative values contract the interval::

    Before expand(-1, -1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                           |<----------------->|
                          end5                end3

    After expand(-1, -1):
    DIS Coordinates:       0   1   2   3   4   5   6   7
                               |<--------->|
                              end5        end3

::

    >>> dis.start, dis.end
    (30, 150)

    >>> # Symmetric expansion
    >>> dis.expand(5).start, dis.expand(5).end
    (25, 155)

    >>> # Asymmetric expansion
    >>> dis.expand(5, 10).start, dis.expand(5, 10).end
    (25, 160)

    >>> # Upstream-only expansion
    >>> dis.expand(5, 0).start, dis.expand(5, 0).end
    (25, 150)

    >>> # Contraction with negative values
    >>> dis.expand(-10, -20).start, dis.expand(-10, -20).end
    (40, 130)

.. note::

    Contracting to exactly zero length is valid, but contracting past
    zero raises ``ValueError``.

Positional Comparisons
======================

``upstream_of`` and ``dnstream_of`` compare two DIS intervals that share
the same coordinate space and the same ``on_coordinate_strand``. Both
methods require strict separation — any overlap returns ``False``.

upstream_of
~~~~~~~~~~~

``upstream_of(other)`` returns ``True`` if ``self`` is strictly 5' of
``other`` with no overlap. Adjacent intervals (where ``self.end`` equals
``other.start``) count as upstream::

    DIS Coordinates:       0   1   2   3   4   5   6   7   8   9
                           |<->|               |<->|
                             a                   b
    a.upstream_of(b) is True    (no overlap)

    DIS Coordinates:       0   1   2   3   4   5   6   7   8   9
                           |<----->|
                              a    |<----->|
                                      b
    a.upstream_of(b) is True    (adjacent: a.end == b.start)

    DIS Coordinates:       0   1   2   3   4   5   6   7   8   9
                           |<--------->|
                              a    |<----->|
                                      b
    a.upstream_of(b) is False   (overlap)

::

    >>> a = DisjointIntervalSequence(coord_ivs, start=10, end=30)
    >>> b = DisjointIntervalSequence(coord_ivs, start=50, end=80)
    >>> a.upstream_of(b)
    True
    >>> b.upstream_of(a)
    False

    >>> # Adjacent intervals count as upstream
    >>> a2 = DisjointIntervalSequence(coord_ivs, start=10, end=50)
    >>> a2.upstream_of(b)
    True

.. note::

    Both intervals must share the same ``coordinate_intervals`` and the
    same ``on_coordinate_strand``, otherwise ``ValueError`` is raised.
    Two zero-length intervals at the same position are neither upstream
    nor downstream of each other.

dnstream_of
~~~~~~~~~~~

``dnstream_of(other)`` is the mirror of ``upstream_of``: it returns
``True`` if ``self`` is strictly 3' of ``other`` with no overlap.
Adjacent intervals count as downstream. The same requirements on shared
coordinate space and strand apply::

    >>> a = DisjointIntervalSequence(coord_ivs, start=50, end=80)
    >>> b = DisjointIntervalSequence(coord_ivs, start=10, end=30)
    >>> a.dnstream_of(b)
    True
    >>> b.dnstream_of(a)
    False

within
~~~~~~

``within(other)`` returns ``True`` if ``self``'s interval is fully
contained within ``other``'s interval. Boundary-inclusive: an interval
is within another if it shares the same start and/or end. An interval
is always within itself. The same requirements on shared coordinate
space and strand apply::

    DIS Coordinates:       0   1   2   3   4   5   6   7   8   9
                                   |<->|
                                     a
                           |<------------->|
                                  b
    a.within(b) is True

    DIS Coordinates:       0   1   2   3   4   5   6   7   8   9
                           |<------------->|
                                  a
                                   |<->|
                                     b
    a.within(b) is False

::

    >>> a = DisjointIntervalSequence(coord_ivs, start=30, end=50)
    >>> b = DisjointIntervalSequence(coord_ivs, start=10, end=80)
    >>> a.within(b)
    True
    >>> b.within(a)
    False

    >>> # An interval is within itself
    >>> a.within(a)
    True

    >>> # Zero-length intervals are within any enclosing interval
    >>> z = DisjointIntervalSequence(coord_ivs, start=50, end=50)
    >>> z.within(a)
    True
