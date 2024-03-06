# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from . import interval as _interval
from ._apply_variants import apply_variants
from ._apply_variants import check_variants_list
from . import _util


class VariantGenome(object):
    """A reference genome with a number of variants applied to it.

    A variant genome behaves like a :py:class:`~genome_kit.Genome` object,
    object, except any requests for DNA sequence (or quantities that depend
    on DNA sequence) will be returned with the variants already applied.
    If the variants change the length of the sequence, then the
    interval's `anchor` and `anchor_offset` determine how the query interval
    will be aligned to the variant genome's coordinate system.

    TODO: see explanation of anchors.
    """

    __slots__ = ('genome', 'variants')

    def __init__(self, reference_genome, variants):
        """Initialize a variant genome.

        Parameters
        ----------
        reference_genome : :py:class:`~genome_kit.Genome`
            The reference genome to which variants are applied.

        variants : :py:class:`~genome_kit.Variant` | :py:class:`list` of (:py:class:`~genome_kit.Variant`)
            Either a single variant or a list of variants. Aligned-normalization of
            variants is not yet implemented, so the user is currently
            responsible for not passing variants that overlap. The variants must be provided
            as a :py:class:`~genome_kit.Variant` object. String variants are no longer allowed,
            please see :py:class:`~genome_kit.Variant` for more information. Since we assume that
            variants never overlap one another this means that all variant
            positions are with respect to the reference genome.
        """

        self.genome = reference_genome
        self.variants = check_variants_list(self.genome, variants)

    def dna(self, interval):
        """Extract variant DNA for an anchored interval or coordinate.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            If `interval` is on the negative strand, the reverse complement
            sequence is returned. The `interval.anchor` and `interval.anchor_offset`
            attributes determine the policy of how length-changing variants are applied.

        Returns
        -------
        :py:class:`str`
            The variant DNA sequence.
        """

        return apply_variants(self.genome.dna, self.variants, interval)

    def find_motif(self, interval, motif, match_position=0, find_overlapping_motifs=False):
        """Find a genomic motif in an interval on a :py:class:`~genome_kit.VariantGenome`.

        This method finds all occurrences of a genomic motif on a VariantGenome. It can
        handle all kinds of variants including indels and complex variants. When a motif
        is found within an insertion, it will set the ``anchor_offset`` attribute on
        the returned :py:class:`~genome_kit.Interval`.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            A genomic interval to search for the motif.
        motif : :py:class:`str`
            A genomic motif to search for.
        match_position : :py:class:`int` | :py:class:`str`
            Very often you want to find a position after a motif or in the middle
            of a motif. For example, when you are searching for a ``AG`` core
            acceptor site motif you want to match the 3p end of the ``AG`` motif to
            return the position of the putative splice site.
            To do this you can specify a match position as string with the values
            ``"3p"`` or ``"5p"``, or you can pass an integer that indicates a position
            in the motif with ``0 <= match_position <= len(motif)``.
        find_overlapping_motifs : :py:class:`bool`
            Whether to return matches for overlapping motifs. The default is to only
            find non-overlapping motifs.

        Returns
        -------
        :py:class`list` of :py:class:`~genome_kit.Interval`
            A list of all matches (length-0 intervals) found in the query interval.

        Example
        -------::
            >>> from genome_kit import Genome, Variant
            >>> genome = Genome('hg19')
            >>> variant = Variant.from_string("chr1:40033:A:G", genome)
            >>> variant_genome = VariantGenome(genome, variant)
            >>> interval = genome.interval('chr1', '+', 40000, 40080)
            >>> motif = 'TAG'
            >>> variant_genome.find_motif(interval, motif)
            [Interval("chr1", "+", 40030, 40030, "hg19", 40030)]
        """
        # TODO: Longterm we'll want to move towards a generator function to be more efficient
        # TODO: on long intervals.

        if match_position == '5p':
            match_position = 0
        elif match_position == '3p':
            match_position = len(motif)

        if not (0 <= match_position <= len(motif)):
            raise ValueError("Require 0 <= match_position <= len(motif) [match_position={}]".format(match_position))

        strand = interval.strand == '+'

        # All motif finding is done on the reference strand. Therefore the motif
        # is reverse complemented when the interval is on the reverse strand. As
        # well, the match_position is recomputed to be applied from the end of the
        # motif (such at it is correct on the reverse strand), e.g. '5p' means
        # match_position = 0 on the forward strand and match_position = len(motif) on
        # the reverse strand.
        if not strand:
            match_position = len(motif) - match_position

        motif = motif.upper()
        dna_start = interval.start

        if not strand:
            motif = _util.reverse_complement(motif)

        # We extend the interval by one to to handle the case when
        # match_position=len(motif) and the motif occurs at the very
        # end of the interval. In that case we would get
        # pos+match_position == len(reference_alignment) which results in
        # an IndexError.
        forward_interval = _interval.Interval(interval.chromosome, '+', interval.start, interval.end + 1,
                                              interval.reference_genome, interval.anchor, interval.anchor_offset)

        # We get both the DNA itself and the reference alignment. The reference alignment
        # is necessary for constructing a Coordinate that is in terms of the reference
        # genome.
        dna, reference_alignment = apply_variants(
            self.genome.dna, self.variants, forward_interval, reference_alignment=True)

        # We truncate only the DNA by one so we don't get any matches past the interval
        dna = dna[:-1]

        motif_hits = []
        offset = 0

        if find_overlapping_motifs:
            offset_step = 1
        else:
            offset_step = len(motif)

        while True:
            pos = dna.find(motif, offset)
            if pos == -1:
                break

            pos_ref = reference_alignment[pos + match_position]

            anchor_offset = 0
            if isinstance(pos_ref, tuple):
                pos_ref, anchor_offset = pos_ref
            elif pos + match_position > 0:
                # if this is immediately after an insertion, we need to set
                # the anchor offset past the insertion so expanding upstream
                # still works
                prev_ref = reference_alignment[pos + match_position - 1]
                if isinstance(prev_ref, tuple):
                    pos_ref, anchor_offset = prev_ref
                    anchor_offset += 1

            anchor = dna_start + pos_ref

            hit = _interval.Interval(
                interval.chromosome,
                interval.strand,
                anchor,
                anchor,
                interval.refg,
                anchor=anchor,
                anchor_offset=anchor_offset)
            motif_hits.append(hit)
            offset = pos + offset_step

        return motif_hits

    def __repr__(self):
        return '<VariantGenome {} with {} variants>'.format(self.reference_genome, len(self.variants))

    def __getattr__(self, name):
        # Intercept any attribute requests that weren't found on the VariantGenome object itself,
        # and forward those requests to the Genome object we're wrapping.

        # TODO: This is problematic: See DGENGINE-1401

        return self.genome.__getattribute__(name)
