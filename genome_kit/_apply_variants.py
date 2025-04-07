# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from . import interval as _interval
from ._util import reverse_complement
from .variant import Variant
from operator import attrgetter
import re

_not_actgn = re.compile(r"[^ACGTN]")


def check_variants_list(reference_genome, variants):
    """Basic input validation for a "list of variants" type argument.

    Checks that each item is a :py:class:`~genome_kit.Variant`,
    and is defined on the correct reference genome.
    Converts a single variant into a list with that variant.
    Does not compare variant `ref` attributes to `reference_genome` sequence.

    Parameters
    ----------
    reference_genome : :py:class:`~genome_kit.Genome`
        The genome on which the variants are expected to be defined.

    variants : :py:class:`~genome_kit.Variant` | :py:class:`list` of (:py:class:`~genome_kit.Variant`)
        Either a single variant or a list of variants.

    Returns
    -------
    :py:class:`list` of :py:class:`~genome_kit.Variant`
        The variants as a list.
    """

    if isinstance(variants, Variant):
        variants = [variants]
    elif not isinstance(variants, (tuple, list)):
        raise TypeError("The 'variants' argument must be a Variant object "
                        "or a list/tuple of Variant objects. Received: {}".format(type(variants)))

    assert (isinstance(variants, (list, tuple)))

    if any(not isinstance(v, Variant) for v in variants):
        raise TypeError("An element of variants is not a Variant object. Got: {}".format(map(type, variants)))

    refg = reference_genome.refg
    for v in variants:
        if v.reference_genome != refg:
            raise ValueError("{}'s genome doesn't match the reference genome.".format(v))
        if _not_actgn.search(v.alt):
            raise ValueError("Invalid ALT sequence '{}'".format(v.alt))

    return variants

def _reverse_tuples_in_alignment(alignment):
    result = []
    tmp = []

    for i in alignment:
        if isinstance(i, tuple):  # If the item is a tuple, add it to the temp list
            tmp.append(i)
        else:
            if len(tmp) > 0:
                # If we encountered a non-tuple and there was a tuple sequence, reverse it
                result.extend(tmp[::-1])
                tmp = []
            result.append(i)

    # If there were any tuples left in temp at the end, reverse and add them
    if len(tmp) > 0:
        result.extend(tmp[::-1])

    return result

def _flip_alignment_for_negative_strand(interval, alignment):
    """Given an alignment for the positive strand, flip it so that it represents the alignment for the negative strand.

    This encompasses two steps: firstly, it flips the provided alignment, since
    negative strand intervals are defined from the 3' to the 5' end. Secondly,
    it translates the alignment so that it is reported in ascending order
    instead of descending order. Here's an example:

    Say we had Interval('chr1', '-', 0, 10, "hg37"), with 2 nt insertion at
    position 5 and a deletion of position 7. If these variants are applied to
    the positive strand, we'd get the following alignment back:

    [0, 1, 2, 3, 4, (5, 0), (5, 1), 5, 6, 8]

    In order to convert this to the negative strand, we'd follow these steps:

    1. Flip the alignment: [8, 6, 5, (5, 1), (5, 0), 4, 3, 2, 1, 0]
    2. Reverse the insertion alignment: [8, 6, 5, (5, 0), (5, 1), 4, 3, 2, 1, 0]
    3. Translate the alignment to be in ascending order:
        3a. [9-8, 9-6, 9-5, (9-5, 0), (9-5, 1), 9-4, 9-3, 9-2, 9-1, 9-0]
        3b. [1, 3, 4, (4, 0), (4, 1), 5, 6, 7, 8, 9]
    4. Increment the tuples by 1 to anchor them to the nucleotide following them:
        [1, 3, 4, (5, 0), (5, 1), 5, 6, 7, 8, 9]
    """
    highest_index = len(interval) - 1
    flipped_alignment = []
    alignment = _reverse_tuples_in_alignment(alignment)
    for i in alignment[::-1]:
        if isinstance(i, tuple):
            flipped_alignment.append((highest_index-i[0]+1, i[1]))
        else:
            flipped_alignment.append(highest_index-i)
    return flipped_alignment


def apply_variants(sequence, variants, interval, reference_alignment=False):
    """Apply variants to a sequence interval.

    This function ignores the interval's strand attribute and always returns
    the sequence on the reference strand.

    Parameters
    ----------
    sequence : :py:class:`~genome_kit.GenomeDNA`
        The source from which to extract DNA.

    variants : :py:class:`list` of :py:class:`~genome_kit.Variant`
        A list of variants.

    interval : :py:class:`~genome_kit.Interval` | :py:class:`~genome_kit.Coordinate`
        An anchored interval/coordinate specifying what to extract.

    reference_alignment : :py:class:`bool`
        Whether to compute the alignment between the variant and reference sequence.

    Returns
    -------
    :py:class:`str`
        The variant DNA sequence.

    :py:class:`list`
        The alignment of the variant sequence with the reference sequence
        (if ``reference_alignment=True``). The alignment is comprised of 
        int/tuple[int, int] offsets to `interval.5p`.
        The tuples denote insertions, where tuple[0] is the offset after the 
        insertion and tuple[1] is the indexes the inserted sequence.

    """

    # TODO: How are anchor_offsets greater than the length of an insertion handled?
    # The anchor offset should just implicitly be truncated to the length of the insertion.
    # Need to test this.

    variants = sorted((vv for v in (x._normalized_variant for x in variants if x.chrom == interval.chrom) for vv in v),
                      key=attrgetter('start'))
    start, end = interval.as_dna0()

    anchor = interval.anchor

    if anchor is None:
        # When anchor is None, the sequence is allowed to change length
        var_sequence = _apply_variants_no_anchor(sequence, variants, interval, reference_alignment)
        if reference_alignment:
            var_sequence, alignment = var_sequence
    elif start < anchor < end:
        # When the anchor falls inside the interval, we split the interval up
        # and apply the variant separately on the left and the right half, then
        # join the sequences together again.

        interval_left = _interval.Interval(interval.chromosome, interval.strand, start, anchor,
                                           interval.reference_genome, interval.anchor, interval.anchor_offset)

        interval_right = _interval.Interval(interval.chromosome, interval.strand, anchor, end,
                                            interval.reference_genome, interval.anchor, interval.anchor_offset)

        var_sequence_left = _apply_variants_right_anchor(sequence, variants, interval_left, reference_alignment)
        var_sequence_right = _apply_variants_left_anchor(sequence, variants, interval_right, reference_alignment)

        if reference_alignment:
            var_sequence_left, alignment_left = var_sequence_left
            var_sequence_right, alignment_right_tmp = var_sequence_right

            # We need to join the reference alignments for the left and the right side,
            # but they both start at zero. Therefore we need to figure out what index
            # we need to add to all the elements in the right reference alignment.
            alignment_right = []

            offset = alignment_left[-1]
            insert_offset = 0
            if isinstance(offset, tuple):
                offset, insert_offset = offset
                insert_offset += 1
            else:
                # handle deletions spanning an anchor
                offset += 1 - (min(0, alignment_left[0]))

            for i in alignment_right_tmp:
                if isinstance(i, tuple):
                    i = (i[0] + offset, i[1] + insert_offset)
                else:
                    i += offset
                    insert_offset = 0
                alignment_right.append(i)
            alignment = alignment_left + alignment_right
        var_sequence = var_sequence_left + var_sequence_right

    elif anchor == start:
        var_sequence = _apply_variants_left_anchor(sequence, variants, interval, reference_alignment)
        if reference_alignment:
            var_sequence, alignment = var_sequence
    elif anchor == end:
        var_sequence = _apply_variants_right_anchor(sequence, variants, interval, reference_alignment)
        if reference_alignment:
            var_sequence, alignment = var_sequence
    elif anchor > end:
        # When the anchor falls outside the interval we just extend the interval up
        # to the anchor and crop it later. TODO: Eventually we'll want to improve this.

        # This method is inefficient since it takes extracts all the sequence
        # up to the anchor.
        start, end = interval.as_dna0()

        tmp_interval = _interval.Interval(interval.chromosome, interval.strand, start, anchor,
                                          interval.reference_genome, interval.anchor, interval.anchor_offset)

        var_sequence = _apply_variants_right_anchor(sequence, variants, tmp_interval, reference_alignment)
        if reference_alignment:
            var_sequence, alignment = var_sequence
            alignment = alignment[:end - start]

        var_sequence = var_sequence[:end - start]
    elif anchor < start:
        start, end = interval.as_dna0()

        tmp_interval = _interval.Interval(interval.chromosome, interval.strand, anchor, end, interval.reference_genome,
                                          interval.anchor, interval.anchor_offset)

        var_sequence = _apply_variants_left_anchor(sequence, variants, tmp_interval, reference_alignment)
        if reference_alignment:
            # When we truncate the reference alignment we need to modify the
            # indices so they start at zero.
            var_sequence, alignment_tmp = var_sequence
            alignment_tmp = alignment_tmp[-(end - start):]

            alignment = []
            for i in alignment_tmp:
                if isinstance(i, tuple):
                    i = (i[0] - alignment_tmp[0], i[1])
                else:
                    i -= alignment_tmp[0]
                alignment.append(i)

        var_sequence = var_sequence[-(end - start):]
    else:
        raise AssertionError("Should never be reached")  # pragma: no cover

    if not interval.is_positive_strand():
        var_sequence = reverse_complement(var_sequence)

        if reference_alignment:
            alignment = _flip_alignment_for_negative_strand(interval, alignment)

    if reference_alignment:
        return var_sequence, alignment
    return var_sequence


def _apply_variants_no_anchor(dna, variants, interval, reference_alignment=False):
    start, end = interval.as_dna0()

    v_interval = _interval.Interval(interval.chrom, '+', start, end, interval.refg)
    variant_sequence = dna(v_interval)

    if reference_alignment:
        alignment = list(range(end - start))

    # Apply the variants from end to avoid reindexing (changing relative start)
    for v in reversed(variants):
        v_start, ref, alt = v.start, v.ref, v.alt
        v_start_rel = v_start - start
        v_end_rel = v_start_rel + len(ref)

        if v_start_rel >= len(variant_sequence):
            continue
        if v_end_rel < 0:
            break

        if v_start_rel < 0:
            # Clip left side of variant
            v_offset = start - v_start
            ref = ref[v_offset:]
            alt = alt[v_offset:]
            v_start_rel = 0

        if v_end_rel > len(variant_sequence):
            # Clip right side of variant
            v_offset = len(variant_sequence) - v_end_rel
            ref = ref[:v_offset]
            alt = alt[:v_offset]
            v_end_rel = len(variant_sequence)

        assert v_end_rel >= 0
        assert v_start_rel >= 0
        assert variant_sequence[v_start_rel:v_end_rel].upper() == ref.upper()
        variant_sequence = variant_sequence[:v_start_rel] + alt + variant_sequence[v_end_rel:]

        if reference_alignment:
            if len(ref) > len(alt):
                # Deletion
                del alignment[v_start_rel + len(alt):v_end_rel]
            elif len(ref) < len(alt):
                # Insertion
                alignment[v_end_rel:v_end_rel] = \
                    [(alignment[v_end_rel], i)
                     for i in range(len(alt) - len(ref))]

    if reference_alignment:
        return variant_sequence, alignment

    return variant_sequence


def _apply_variants_right_anchor(dna, variants, interval, reference_alignment=False):
    start, end = interval.as_dna0()
    original_start = start

    variants_processed = []
    # Search from anchor to start for variants extended before original interval
    for v in reversed(variants):
        v_start, ref, alt = v.start, v.ref, v.alt
        v_end = v_start + len(ref)

        # Variant falls outside of interval
        if v_start >= end and v_end > end:
            continue
        if v_end <= start:
            break

        if end == v_start and not ref:
            alt = alt[:interval.anchor_offset]

        # Clip variant if it overlaps the anchor
        if v_start < end <= v_end:
            v_len = end - v_start
            ref = ref[:v_len]
            alt = alt[:v_len]

        # Clip if variant overlaps the start
        if v_start - len(alt) + len(ref) < start:
            v_offset = start - (v_start - len(alt) + len(ref))
            v_start += min(v_offset, len(ref))
            ref = ref[v_offset:]
            alt = alt[v_offset:]

        # Move the start to accommodate the variant
        start += len(alt) - len(ref)
        variants_processed.append((v_start, ref, alt))

    if reference_alignment:
        # need an extra element for insertion at end
        alignment = list(range(start - original_start, end - min(start, original_start) + 1))

    v_interval = _interval.Interval(interval.chrom, '+', start, end, interval.refg)
    variant_sequence = dna(v_interval)

    # Apply the variants from end to avoid reindexing (changing relative start)
    for v_start, ref, alt in variants_processed:
        v_start_rel = v_start - start
        v_end_rel = v_start_rel + len(ref)

        assert v_end_rel >= 0
        assert v_start_rel >= 0
        assert variant_sequence[v_start_rel:v_end_rel].upper() == ref.upper()
        variant_sequence = variant_sequence[:v_start_rel] + alt + variant_sequence[v_end_rel:]

        if reference_alignment:
            if len(ref) > len(alt):
                # Deletion
                del alignment[v_start_rel + len(alt):v_end_rel]
            elif len(ref) < len(alt):
                # Insertion
                alignment[v_end_rel:v_end_rel] = \
                    [(alignment[v_end_rel], i)
                     for i in range(len(alt) - len(ref))]

    variant_sequence = variant_sequence[-len(interval) :]
    assert len(variant_sequence) == len(interval)

    if reference_alignment:
        alignment = alignment[:len(variant_sequence)]
        return variant_sequence, alignment

    return variant_sequence


def _apply_variants_left_anchor(dna, variants, interval, reference_alignment=False):
    start, end = interval.as_dna0()
    original_end = end

    variants_processed = []
    # Search from anchor to end  for variants extended after original interval
    for v in variants:
        v_start, ref, alt = v.start, v.ref, v.alt
        if start == v_start and not ref:
            # Handle anchor_offset cases separately
            alt = alt[interval.anchor_offset:]
            # Clip if variant overlaps the end
            if v_start + len(alt) - len(ref) > end:
                v_clip = v_start + len(alt) - len(ref) - end
                ref = ref[:-v_clip]
                alt = alt[:-v_clip]
            variants_processed.append((v_start, ref, alt))
            continue

        v_end = v_start + len(ref)

        # Variant falls outside of interval
        if v_end <= start:
            continue
        if v_start >= end:
            break

        if v_start <= start < v_end:
            # Variant overlaps the anchor
            v_offset = start - v_start
            ref = ref[v_offset:]
            alt = alt[v_offset:]
            v_start = start

        # Clip if variant overlaps the end
        if v_start + len(alt) - len(ref) > end:
            v_clip = v_start + len(alt) - len(ref) - end
            ref = ref[:-v_clip]
            alt = alt[:-v_clip]

        variants_processed.append((v_start, ref, alt))
        end -= len(alt) - len(ref)

    if reference_alignment:
        alignment = list(range(max(original_end, end) - start))

    v_interval = _interval.Interval(interval.chrom, '+', start, end, interval.refg)
    variant_sequence = dna(v_interval)

    # Apply the variants from end to avoid reindexing (changing relative start)
    for v_start, ref, alt in reversed(variants_processed):
        v_start_rel = v_start - start
        v_end_rel = v_start_rel + len(ref)

        assert v_start_rel <= len(variant_sequence)

        if v_end_rel > len(variant_sequence):
            # Clip variant if it extends beyond the sequence
            v_clip = len(variant_sequence) - v_start_rel
            ref = ref[:v_clip]
            alt = alt[:v_clip]
            v_end_rel = len(variant_sequence)

        assert variant_sequence[v_start_rel:v_end_rel].upper() == ref.upper()
        variant_sequence = variant_sequence[:v_start_rel] + alt + variant_sequence[v_end_rel:]

        if reference_alignment:
            if len(ref) > len(alt):
                # Deletion
                del alignment[v_start_rel + len(alt):v_end_rel]
            elif len(ref) < len(alt):
                # Insertion
                alignment[v_end_rel:v_end_rel] = \
                    [(alignment[v_end_rel], i)
                     for i in range(len(alt) - len(ref))]

    variant_sequence = variant_sequence[: len(interval)]
    assert len(variant_sequence) == len(interval)

    if reference_alignment:
        alignment = alignment[: len(interval)]
        return variant_sequence, alignment

    return variant_sequence
