.. _diseq:

-------------------------------
Disjoint Interval Sequences
-------------------------------

Motivation
==========

When working with transcripts, the situation may arise where we want to ignore the
introns (or other features) of the transcript. If you were to represent the transcript
with those parts removed, you would be left with a disjoint series of :py:class:`~genome_kit.Interval`
objects that are difficult to work with directly (for example, creating an interval on
this disjoint space, or querying the position of a specific sequence within a
CDS, relative to the spliced RNA sequence).

For this reason, the :py:class:`~genome_kit.diseq.DisjointIntervalSequence` class was
introduced. :py:class:`~genome_kit.diseq.DisjointIntervalSequence` (DIS) simplifies
working with intervals that exist on a disjoint coordinate space.


Overview
========

A :py:class:`~genome_kit.diseq.DisjointIntervalSequence` (DIS) represents
a flattened coordinate system over a sequence of disjoint genomic intervals.
For example, the exons of a transcript form a disjoint set of genomic intervals that,
when concatenated, represent the spliced RNA sequence.

A DIS has two layers:

- A **coordinate space**: the underlying genomic
  :py:class:`~genome_kit.Interval` objects (e.g. exons) that define the
  flattened index system. These intervals are sorted 5'→3' and must not overlap.
- An **interval**: a sub-range within that coordinate space, defined by
  a start and end index, where start <= end.

To explain how the coordinate space and interval layers interact, let's ignore code for
now, and just use some diagrams to illustrate the concepts.

Say we have a transcript on the + strand represented by the diagram below:
::
    Genomic Coordinates:  153 154 155 156 157 158 159 160 161 162 163 164 165 166 167
    DNA Sequence:          A   T   G   C   C   G   C   A   T   G   C   C   G   C
                          |<------->| |<------->| |<--->| |<----------->| |<--->|
                             Exon1      Intron1    Exon2      Intron2      Exon3

If we were to take only the exons, we would have the following disjoint intervals:
::
    Genomic Coordinates:  153 154 155 159 160 165 166 167
    DNA Sequence:          A   T   G   C   A   G   C
                          |<------->| |<--->| |<--->|
                             Exon1     Exon2   Exon3

Let's say we want to create an interval on this series of disjoint exon intervals,
spanning from the start of Exon1 to the end of Exon3. We can start by converting our
list of exons into a DIS coordinate space
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          A   T   G   C   A   G   C
                           |<----->|   |<->|   |<->|
                             Exon1     Exon2   Exon3

Now let's place the interval on the DIS coordinate space
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          A   T   G   C   A   G   C
                           |<--------------------->|
                          end5      Interval      end3
    Start Index:     0
    End Index:       7

We see that the interval spans the full length of the coordinate space, and is defined by
a start index of 0 and an end index of 7.

.. note::
    Why 7 and not 6? The disjoint interval follows the convention of
    :py:class:`~genome_kit.Interval` where intervals are half-open
    (the end index is exclusive).

The above example illustrates the basics of how a DIS works. However, it is possible to
do more. We can instead define an interval within the DIS on the strand opposite that 
of the coordinate space. Let's start with the DIS coordinate space from above
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          A   T   G   C   A   G   C
                           |<----->|   |<->|   |<->|
                             Exon1     Exon2   Exon3

Now let's add the negative (opposite) strand to the diagram
::
                                Positive Strand
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence (+):      A   T   G   C   A   G   C
    -----------------------------------------------------
    DNA Sequence (-):      T   A   C   G   T   C   G
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                Negative Strand

Importantly, we see the DIS coordinates are the same on both strands. This simplifies
things when you want to get the complement of a given interval, as you can use the same
indices and just flip the strand. To illustrate this, let's now define the same interval
as before (spanning the entire coordinate space) but on the negative strand
::
                                Positive Strand
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence (+):      A   T   G   C   A   G   C
    -----------------------------------------------------
    DNA Sequence (-):      T   A   C   G   T   C   G
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                Negative Strand
                           |<--------------------->|
                          end3      Interval      end5
    Start Index:     0
    End Index:       7
    On Coordinate Strand: False

You will notice "On Coordinate Strand: False" has been added to the diagram. Since we
aren't able to determine which strand the interval is on just from the indices, this
variable is used to let us know the strandedness of the interval.

Thus far we have only defined a DIS from intervals on the + strand. When defining
a DIS from intervals on the negative strand, much remains the same. However, there is
one important difference from a regular Interval: On a DIS created from negative-strand
intervals, the indices still increase in the 5'→3' direction of the transcript. Let's
take a look at an example:

Say we have the following transcript on the negative strand represented by the diagram
below:
::
    Reminder: On the - strand, the 3' end is on the left and the 5' end is on the right!

                             Exon3      Intron2    Exon2      Intron1      Exon1
                          |<------->| |<------->| |<--->| |<----------->| |<--->|
    DNA Sequence (-):      G   T   C   A   G   T   C   A   G   T   C   A   G   T
    Genomic Coordinates:  153 154 155 156 157 158 159 160 161 162 163 164 165 166 167
                                            Negative Strand
Taking just the exons:
::
                             Exon3     Exon2   Exon1
                          |<------->| |<--->| |<--->|
    DNA Sequence (-):      G   T   C   C   A   G   T
    Genomic Coordinates:  153 154 155 159 160 165 166 167
                                    Negative Strand

Now let's create a DIS from these exons:
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
                           |<->|   |<->|   |<----->|
                           Exon1   Exon2     Exon3

Notice that the sequence flips relative to the direction of the indices. What has
happened is that the DIS coordinate space is defined in the 5'→3' direction of the
transcript, regardless of genomic strand. In a DIS, 0 always corresponds to the
DIS coordinate's 5' end, and the largest index corresponds to the DIS coordinate's 3'
end.

Let's now define an interval on this DIS
::
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
                           |<--------------------->|
                          end5      Interval      end3
    Start Index:     0
    End Index:       7
    On Coordinate Strand: True

We see that despite creating the DIS from the negative strand, the full-length interval
on the coordinate strand still looks the same as in the + strand example. When working
with DIS objects, you only need to think of things in terms of "same strand" or
"opposite strand".

To complete the example, let's define an interval on this DIS that is on the opposite strand of the coordinate space
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

Now that you understand how a DIS works conceptually, you can read on to see how to
manipulate them in code.

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
You can also specify a region to use CDS or UTR intervals::

    >>> dis_cds = DisjointIntervalSequence.from_transcript(transcript, region="cds")
    >>> dis_utr5 = DisjointIntervalSequence.from_transcript(transcript, region="utr5")
    >>> dis_utr3 = DisjointIntervalSequence.from_transcript(transcript, region="utr3")

The ``coord_id`` and ``interval_id`` default to ``transcript.id`` but can
be overridden::

    >>> dis = DisjointIntervalSequence.from_transcript(
    ...     transcript, coord_id="my_coord", interval_id="my_interval")

From Intervals
~~~~~~~~~~~~~~

You can also construct a DIS from any sequence of
:py:class:`~genome_kit.Interval` objects (or annotation objects like
:py:class:`~genome_kit.Exon` that have an ``.interval`` attribute)::

    >>> from genome_kit import Interval
    >>> exon_intervals = [e.interval for e in transcript.exons]
    >>> dis = DisjointIntervalSequence.from_intervals(exon_intervals, coord_id="my_coord")

The intervals must all share the same chromosome, strand, and reference
genome. They are automatically sorted 5'→3' and checked for overlaps.

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
    >>> dis.coord_transcript_strand
    '+'
    >>> dis.reference_genome
    'hg38'
    >>> dis.coord_id
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
interval's strand::
    On coordinate strand (on_coordinate_strand=True):
    Start Index:     1
    End Index:       6
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
                               |<------------->|
    -----------------------------------------------------
    DNA Sequence:          A   C   T   G   G   A   C
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                  Opposite Strand

    >>> dis = DisjointIntervalSequence.from_transcript(transcript)
    >>> dis.end5_index   # same as start when on coordinate strand
    1
    >>> dis.end3_index   # same as end when on coordinate strand
    6


    Off coordinate strand (on_coordinate_strand=False):
    Start Index:     3
    End Index:       7
    DIS Coordinates:       0   1   2   3   4   5   6   7
    DNA Sequence:          T   G   A   C   C   T   G
    -----------------------------------------------------
    DNA Sequence:          A   C   T   G   G   A   C
    DIS Coordinates:       0   1   2   3   4   5   6   7
                                       |<--------->|
                                  Opposite Strand

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
