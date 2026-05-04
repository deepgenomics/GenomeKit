.. _anchors:

Anchors
=======

When applying variants to a genome and extracting features (most importantly DNA sequence) from the resulting
:py:class:`~genome_kit.VariantGenome`, there needs to be some notion of alignment between the wild-type and mutant
sequences when variants change the length of the sequence.  Furthermore, many machine learning methods expect
fixed-length inputs so there needs to a predictable and unambiguous mechanism of how variants affect the alignment.

GenomeKit solves this problem through a mechanism called *anchors*. Rather than allowing arbitrary alignments, the
alignment is fixed at a single coordinate. Usually this is a position of interest such as a splice site. It is
important to understand the anchor concept to understand the objects and sequences returned by GenomeKit, but in
practice it is rarely necessary to set anchors manually as functions that return
:py:class:`~genome_kit.Interval` typically set the anchor appropriately.


Conceptually, an anchor is a coordinate (a *"gap"* between two bases) that will always align between the wild-type
and mutant sequence. Consider the interval ``Interval("chr7", "+", 117231977, 117231997, "hg19")`` which is a short
window of size 20 around the acceptor site of exon 13 on CFTR::

    117,231,977                   117,231,987                117,231,996
    0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
    T  A  T  C  T  T  A  A  A  G  C  T  G  T  G  T  C  T  G  T

The first line shows the genomic coordinates of the first and last elements in the sequence as well as the coordinate
of the splice site in the middle of the interval.

.. note::
   In this chapter we stick to the GenomeKit convention of zero-based genomic coordinates. However, any coordinates in
   variants strings like ``'chr7:117,231,993:TCT:T'`` are in one-based coordinates as in the Clinvar variant format.

Now consider we are going to apply a length changing variant to this interval. Consider the variant
``chr7:117,231,993:TCT:T`` which deletes the ``CT`` dinucleotide at position 16.


Intervals without anchors
-------------------------

The default behaviour of :py:meth:`genome_kit.VariantGenome.dna` is to apply the variant naïvely (without an anchor)
and allow the length of the sequence to change::

    >>> genome = Genome('hg19')
    >>> variant_genome = VariantGenome(genome, 'chr7:117,231,993:TCT:T')
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19")
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'TATCTTAAAGCTGTGTGT'
    >>> len(var_seq)
    18

As you can see two bases were removed and the resulting mutant sequence is shorter than the wild-type.


Deletions with an anchored interval
-----------------------------------

If we set an anchor, we can specify to keep the wild-type and mutant sequences aligned at the splice site and expand
the mutant sequence to keep its length constant::

    >>> genome = Genome('hg19')
    >>> variant_genome = VariantGenome(genome, 'chr7:117,231,993:TCT:T')
    >>> # Anchor set at the splice site
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", anchor=117231987)
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'TATCTTAAAGCTGTGTGTAA'
    >>> len(var_seq)
    20

What happened here is best understood when looking at the aligned wild-type and mutant sequences::

        0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19
    WT: T  A  T  C  T  T  A  A  A  G | C  T  G  T  G  T  C  T  G  T
    MT: T  A  T  C  T  T  A  A  A  G | C  T  G  T  G  T  .  .  G  T  A  A
        0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 .  .  16 17 18 19

The vertical bars indicate the anchor. Since the deletion occurred on the right side of the anchor, the interval is
extended on the right side to include extra sequence equal to the length of the deletion.

Specifying anchors
------------------

An anchor is always an integer coordinate on the chromosome. The anchor can be at the start or end of an interval,
anywhere inside the interval, or even outside the interval (we'll look at this case in :ref:`insertions-at-the-anchor`).
However, :py:class:`~genome_kit.Interval` allows to set the anchor by passing one of the strings ``'5p'``, ``'3p'``, or
``'center'`` to set the anchor more conveniently.


Insertions with an anchored interval
------------------------------------

Now we consider an insertion to the left side of the anchor. Whereas the deletion *pulled in* extra sequence to keep
the length constant, the insertion will *push out* sequence from the interval::

    >>> genome = Genome('hg19')
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", anchor=117231987)
    >>> # Insert GGG at 117,231,982
    >>> variant_genome = VariantGenome(genome, 'chr7:117,231,982:T:TGGG')
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'CTGGGTAAAGCTGTGTCTGT'
    >>> len(var_seq)
    20

What happens is more easily seen when looking at the alignment::

        0  1  2  3  4  .  .  .  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19
    WT: T  A  T  C  T  .  .  .  T  A  A  A  G | C  T  G  T  G  T  C  T  G  T
    MT: .  .  .  C  T  G  G  G  T  A  A  A  G | C  T  G  T  G  T  C  T  G  T
        .  .  .  0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19

As in the first example, the anchor coordinates stay aligned and the left or right side of the interval grow or shrink
to accommodate the length change.


Deletions at the anchor
-----------------------

Since the anchor is an idealized coordinate or *"gap"*, the mechanism still works even when a deletion occurs right
around the anchor. For example, if a five deletion overlaps the anchor and three bases are to the left of the anchor
and two bases are to the right of the anchor, then the interval start will be extended by three bases and the end by
two bases to accommodate the deletion. Coming back to our previous example::

    >>> genome = Genome('hg19')
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", anchor=117231987)
    >>> # Delete five bases at 117,231,984
    >>> variant_genome = VariantGenome(genome, 'chr7:117,231,984:AAAGCT:A')
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'TTATATCTTAGTGTCTGTAA'

Let's examine the concrete alignment again, so see what happened::

                 0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19
    WT:          T  A  T  C  T  T  A  A  A  G | C  T  G  T  G  T  C  T  G  T
    MT: T  T  A  T  A  T  C  T  T  A  .  .  . | .  .  G  T  G  T  C  T  G  T  A  A
        0  1  2  3  4  5  6  7  8  9  .  .  . | .  .  10 11 12 13 14 15 16 17 18 19

As explained above the five deletions were split between the left and right side to the anchor and the sequence was
expanded accordingly to keep the length constant.

.. _insertions-at-the-anchor:

Insertions at the anchor
------------------------

The last case is slightly more complicated. When insertions occur directly at the anchor, it is not clear whether
the insertion should occur on the left or ride side of the anchor. There is also the special case when we would like to
align a coordinate *within the insertion*. This case is actually important to distinguish, *i.e.* when an insertion
contains an entire core splice site motif and we would like to align this de-novo splice site with a coordinate on
the reference genome. To handle of these cases, there is an additional argument to :py:class:`genome_kit.Interval`
and :py:class:`genome_kit.VariantGenome`, called the ``anchor_offset`` that specifies *where to align within an
insertion*.

First consider the default case. If we do not pass the argument, by default ``anchor_offset=0`` and **insertions will
occur right of the anchor**. Here is an example where we insert ``TTT`` directly at the anchor coordinate::

    >>> genome = Genome('hg19')
    >>> # Insert TTT at 117,231,987
    >>> # Implies anchor_offset=0
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", 117231987)
    >>> variant_genome = VariantGenome(genome, 'chr7:117231987:G:GTTT')
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'TATCTTAAAGTTTCTGTGTC'

When you study the alignment you can see that the insertion occurred on the right side of the anchor::

        0  1  2  3  4  5  6  7  8  9 | .  .  .  10 11 12 13 14 15 16 17 18 19
    WT: T  A  T  C  T  T  A  A  A  G | .  .  .  C  T  G  T  G  T  C  T  G  T
    MT: T  A  T  C  T  T  A  A  A  G | T  T  T  C  T  G  T  G  T  C
        0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19

The opposite case is when we want to insert on the left of the anchor instead, which we can achieve by setting
``anchor_offset=len(insertion)``. Using the same variant as before we now get a different result when we set the
anchor offset::

    >>> genome = Genome('hg19')
    >>> # Insert TTT at 117,231,987
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", 117231977, anchor_offset=3)
    >>> variant_genome = VariantGenome(genome, 'chr7:117231987:G:GTTT')
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'CTTAAAGTTTCTGTGTCTGT'

By examining the alignment we can see that the insertion now occurred on the left side of the anchor::

        0  1  2  3  4  5  6  7  8  9  .  .  . | 10 11 12 13 14 15 16 17 18 19
    WT: T  A  T  C  T  T  A  A  A  G  .  .  . | C  T  G  T  G  T  C  T  G  T
    MT: .  .  .  C  T  T  A  A  A  G  T  T  T | C  T  G  T  G  T  C  T  G  T
                 0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19

Of course we can also set the anchor offset to other values. In the above example if we want one base to be
inserted to the left and two bases to be inserted to the right of the anchor we can set ``anchor_offset=1``::

    >>> genome = Genome('hg19')
    >>> # Insert TTT at 117,231,987
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", 117231977, anchor_offset=1)
    >>> variant_genome = VariantGenome(genome, 'chr7:117231987:G:GTTT')
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'ATCTTAAAGTTTCTGTGTCT'

The alignment confirms this is what happened::

        0  1  2  3  4  5  6  7  8  9  . | .  .  10 11 12 13 14 15 16 17 18 19
    WT: T  A  T  C  T  T  A  A  A  G  . | .  .  C  T  G  T  G  T  C  T  G  T
    MT: .  A  T  C  T  T  A  A  A  G  T | T  T  C  T  G  T  G  T  C  T
           0  1  2  3  4  5  6  7  8  9 | 10 11 12 13 14 15 16 17 18 19

.. note::

   In general, anchors are not specific to a particular variant. An interval or coordinate with a set anchor and/or
   anchor offset can be applied to any reference or variant genome. Anchor offsets are most often used to achieve an
   alignment of a coordinate in a particular insertion with the reference genome but are in principle independent of
   any variant, It is more likely that unintended results are returned when intervals that have an anchor offset are
   applied to different variants. Therefore, the following should be kept in mind:

   1. Anchor offsets only have an effect when applying an insertion directly at the anchor. An insertion at any other
      position is not affected by the anchor offset.
   2. The effective anchor offset cannot be greater than the length of the insertion. For example, if an anchor offset
      of ten is set and the interval is given an insertion of length five, then only an anchor offset of five is
      applied, *i.e.* the insertion is done on the left side of the anchor (so you could even set the anchor offset to
      ``MAX_INT`` if you always wanted to insert on the left of the anchor).


.. _anchors-outside-an-interval:

Anchors outside an interval
---------------------------

It is even possible to define an anchor outside of an interval. This is primarily done for completeness to allow an
interval to have its length or coordinates modified, but still work for feature extraction when the anchor stays the
same. Conceptually, defining an anchor outside an interval is equivalent to extending the interval up to the anchor,
extracting the sequence, and then trimming the returned sequence by the same amount it was first extended.

When the anchor falls outside the interval, it is possible that indels outside the interval affect the returned
sequence. Consider our previous interval, but now we set the anchor at 117,231,967, ten bases before the start of the
interval. Then we apply a deletion between the anchor and the start of the interval::

    >>> genome = Genome('hg19')
    >>> # Insert TTT at 117,231,987
    >>> interval = Interval("chr7", "+", 117231977, 117231997, "hg19", anchor=117231967)
    >>> variant_genome = VariantGenome(genome, 'chr7:117231970:TGT:T')
    >>> var_seq = variant_genome.dna(interval)
    >>> var_seq
    'TCTTAAAGCTGTGTCTGTAA'

Again, let's check alignment. Here we denote the anchor by the vertical line and the interval by parentheses::

       |-10 -9 -8 -7 -6 -5 -4 -3 -2 -1[0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19]
    WT:| T  A  T  G  T  T  T  T  T  A [T  A  T  C  T  T  A  A  A  G C  T  G  T  G  T  C  T  G  T ]
    MT:| T  A  T  .  .  T  T  T  T  A  T  A [T  C  T  T  A  A  A  G C  T  G  T  G  T  C  T  G  T  A  A ]
                                            [0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19]

The returned sequence is the same as if we had first extended the start of the interval up to the anchor by ten
bases, then applied the interval to the variant genome and trimmed the first ten bases off the start of the returned
sequence again.
