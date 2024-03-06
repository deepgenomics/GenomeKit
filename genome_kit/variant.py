# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import pickle
from bisect import bisect_left
from builtins import zip
from operator import attrgetter

from . import _cxx
from ._cxx_util import mock, mock_result, mock_unreachable, strip_mock_bases
from .interval import Interval

_NUCLEOTIDE_SYMBOLS = 'ACGTWSMKRYBDHVNZ'


@strip_mock_bases
@_cxx.register
class Variant(_cxx.Variant, Interval):
    """A variant.

    Bases: :py:class:`~genome_kit.Interval`
    """
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    # noinspection PyMissingConstructor
    @mock
    def __init__(self, chromosome, start, ref, alt, reference_genome):  # pragma: no cover
        """Initialize a variant from a DNA0 position.

        Parameters
        ----------
        chromosome : :py:class:`str`
            The chromosome name ('chr1', ...).

        start : :py:class:`int`
            The start of the reference genome spanned by `ref` (inclusive, 0-based).

        ref : :py:class:`str`
            The reference allele.

        alt : :py:class:`str`
            The alternate allele.

        reference_genome : :py:class:`str` | :py:class:`~genome_kit.Genome`
            The reference genome. See :py:class:`~genome_kit.Interval`.
        """
        mock_unreachable()

    @mock
    @property
    def interval(self):  # pragma: no cover
        """The interval spanned by this variant in the reference genome.

        Note that :py:class:`~genome_kit.Variant` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this variant's `ref` sequence in the reference genome.
        """
        return mock_result(Interval)

    @mock
    @property
    def position(self):  # pragma: no cover
        """The integer position of the start (inclusive, 0-based).

        This property is a synonym for `start` and is provided for
        historical reasons.  For consistency with other interval-based
        objects, prefer the `start` property.

        Returns
        -------
        :py:class:`int`
            The start position of this interval within the reference genome.

        """
        return mock_result(int)

    @mock
    @property
    def ref(self):  # pragma: no cover
        """The reference allele as a DNA string.

        Returns
        -------
        :py:class:`str`
            The sequence of the reference allele that is altered by this variant.
        """
        return mock_result(str)

    @mock
    @property
    def alt(self):  # pragma: no cover
        """The alternate allele as a DNA string.

        Returns
        -------
        :py:class:`str`
            The sequence of the alternate allele that this variant represents.
        """
        return mock_result(str)

    def as_variant_string(self):
        """Returns a string representation of the variant.

        See Also
        --------
        :py:meth:`~genome_kit.Variant.from_string`

        Returns
        -------
        :py:class:`str`
            The variant in the format ``"chromosome:start:ref:alt"``, where
            start is in DNA1 coordinates.
        """
        return '{}:{}:{}:{}'.format(self.chromosome, self.start + 1, self.ref, self.alt)

    def __repr__(self):
        return '<Variant {}:{}:{}:{}:{}>'.format(self.chrom, self.start + 1, self.ref, self.alt, self.refg)

    __str__ = as_variant_string

    @staticmethod
    def from_string(variant, reference_genome):
        """Initialize a variant object from a variant string.

        Parameters
        ----------
        variant : :py:class:`str`
            A variant in the format ``"chromosome:start:ref:alt"`` where start is in DNA1 coordinates.

        reference_genome : :py:class:`~genome_kit.Genome`
            The reference genome that `variant` diverges from.

        See Also
        --------
        :py:meth:`~genome_kit.Variant.as_variant_string`

        """
        try:
            components = variant.split(':')
        except (AttributeError, ValueError, SyntaxError):
            raise ValueError("Invalid variant: {}.".format(variant))

        if not len(components) == 4:
            raise ValueError("Variant must have exactly four components: {}".format(variant))

        chromosome, position, ref, alt = components

        try:
            # Allow position to have comma separators
            position = int(position.replace(',', '')) - 1
        except ValueError:
            raise ValueError("Invalid position in variant: {}".format(variant))

        chromosome, ref, alt = Variant._preprocess_variant(chromosome, ref, alt)
        variant_obj = Variant(chromosome, position, ref, alt, reference_genome)
        try:
            variant_obj._validate_variant(reference_genome)
        except ValueError as e:
            raise ValueError(f"Variant string was invalid: {variant}") from e
        return variant_obj

    @staticmethod
    def spanning(x, y, variants=None):
        """Extends :meth:`~genome_kit.Interval.spanning` such that it is aware
        of variants and anchors.

        Example::

            >>> from genome_kit import Genome, Variant, Interval
            >>> hg19 = Genome("hg19")
            >>> v = Variant("chr1", 100, "", 10 * "N", hg19)
            >>> a = Interval("chr1", "+", 100, 150, hg19, 100)
            >>> b = Interval("chr1", "+", 200, 250, hg19)
            >>> Interval.spanning(a, b)
            Interval("chr1", "+", 100, 250, "hg19")
            >>> Variant.spanning(a, b, [v])
            Interval("chr1", "+", 100, 260, "hg19", 100)

        Parameters
        ----------
        x : :py:class:`~genome_kit.Interval`
            an interval

        y : :py:class:`~genome_kit.Interval`
            another interval

        variants : :class:`list` of :class:`~genome_kit.Variant` or None, optional
            variants to apply to `x` and `y`; variants cannot overlap.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The resulting spanning interval.
        """

        if x.anchor is None and y.anchor is None:
            return Interval.spanning(x, y)

        if not variants:
            raise ValueError("Anchored intervals (%s, %s) require a list of variants." % (x, y))

        # verify input
        chrom = x.chrom
        if chrom != y.chrom:
            raise ValueError("Both intervals must be on the same chromosome, got: %s, %s" % (chrom, y.chrom))

        strand = x.strand
        if strand != y.strand:
            raise ValueError("Both intervals must be on the same strand, got: %s, %s" % (strand, y.strand))

        refg = x.refg
        if refg != y.refg:
            raise ValueError(
                "Both intervals must on on the same reference genome, got: %s, %s" % (refg, y.refg))

        if x.anchor is None or (y.anchor is not None and (
                x.anchor > y.anchor or x.anchor == y.anchor and x.anchor_offset > y.anchor_offset)):
            # simplify with x being the leftmost anchor
            x, y = y, x

        anchor = x.anchor
        anchor_offset = x.anchor_offset
        start, end = x.as_dna0()
        if anchor < start or anchor > end:
            raise ValueError("Outside anchors are unsupported: %s" % x)

        def clamp_anchor_offset(variants, starts, anchor, anchor_offset):
            i = bisect_left(starts, anchor)
            if i != len(starts) and starts[i] == anchor:
                v = variants[i]
                if not v.ref:
                    return min(anchor_offset, len(v.alt))
            return 0

        variants = sorted((vv for n in (v._normalized_variant for v in variants if v.chrom == chrom) for vv in n),
                          key=attrgetter('start'))
        starts = None
        if anchor_offset > 0:
            starts = [v.start for v in variants]
            anchor_offset = clamp_anchor_offset(variants, starts, anchor, anchor_offset)

        # map y into x's anchored interval and then union
        if y.anchor is None:
            if y.start < anchor:
                length = Variant._get_variant_length(variants, Interval(chrom, strand, y.start, anchor, refg))
                length += anchor_offset
                start = min(start, anchor - length)
            if y.end > anchor:
                length = Variant._get_variant_length(variants, Interval(chrom, strand, anchor, y.end, refg))
                length -= anchor_offset
                end = max(end, anchor + length)
            return Interval(chrom, strand, start, end, refg, anchor, x.anchor_offset)

        # map y's anchor (such that x.anchor + length = y.anchor)
        if y.anchor < y.start or y.anchor > y.end:
            raise ValueError("Outside anchors are unsupported: %s" % y)

        length = Variant._get_variant_length(variants, Interval(chrom, strand, anchor, y.anchor, refg))
        length -= anchor_offset
        if y.anchor_offset > 0:
            if not starts:
                starts = [v.start for v in variants]
            length += clamp_anchor_offset(variants, starts, y.anchor, y.anchor_offset)

        delta_y_anchor = anchor + length - y.anchor
        start = min(start, delta_y_anchor + y.start)
        end = max(end, delta_y_anchor + y.end)
        return Interval(chrom, strand, start, end, refg, anchor, x.anchor_offset)

    @property
    def _normalized_variant(self):
        ref = self.ref
        alt = self.alt
        overlap = min(len(ref), len(alt))

        # strip padding
        left_pad = next((i for i, (x, y) in enumerate(zip(ref, alt)) if x != y), overlap)
        ref = ref[left_pad:]
        alt = alt[left_pad:]
        overlap -= left_pad
        right_pad = next((i for i, (x, y) in enumerate(zip(reversed(ref), zip(reversed(alt))))), overlap)
        ref = ref[:len(ref) - right_pad]
        alt = alt[:len(alt) - right_pad]

        start = self.start + left_pad
        if len(alt) == len(ref):
            return [] if not ref else [self if ref == self.ref else self._with_ref(start, ref, alt)]

        # Split up complex variants into a substitution
        # plus an insertion or deletion
        # Also used to normalize Clinvar variants
        if len(alt) > len(ref):
            # Split into substitution and insertion
            insertion = self._with_ref(start + len(ref), "", alt[len(ref):])
            return [insertion] if not ref else [self._with_ref(start, ref, alt[:len(ref)]), insertion]

        # Split into substitution and deletion
        deletion = self._with_ref(start + len(alt), ref[len(alt):], "")
        return [deletion] if not alt else [self._with_ref(start, ref[:len(alt)], alt), deletion]

    def _with_ref(self, start, ref, alt):
        return Variant(self.chrom, start, ref, alt, self.refg)

    @staticmethod
    def _preprocess_variant(chromosome, ref, alt):
        if not chromosome.startswith('chr'):
            chromosome = 'chr' + chromosome

        ref = ref.upper()
        alt = alt.upper()

        if ref in '-.':
            ref = ""
        if alt in '-.':
            alt = ""

        return chromosome, ref, alt

    def _validate_variant(self, genome):
        for r in self.ref:
            if r not in _NUCLEOTIDE_SYMBOLS:
                raise ValueError("Invalid nucleotide in ref allele: {}".format(repr(self.ref)))

        for a in self.alt:
            if a not in _NUCLEOTIDE_SYMBOLS:
                raise ValueError("Invalid nucleotide in alt allele: {}".format(repr(self.alt)))

        interval = Interval(self.chromosome, '+', self.position, self.position + len(self.ref), genome)
        wt_seq = genome.dna(interval)

        if not wt_seq == self.ref:
            raise ValueError("Variant's ref sequence does not match the genome ({} != {})".format(wt_seq, self.ref))

    @staticmethod
    def _get_variant_length(sorted_normalized_variants, interval):
        start = interval.start
        length = len(interval)

        # Match ordering used in _apply_variants_no_anchor
        for v in reversed(sorted_normalized_variants):
            v_start_rel = v.start - start
            if v_start_rel >= length:
                continue

            v_end_rel = v_start_rel + len(v.ref)
            if v_end_rel < 0:
                break

            alt_length = len(v.alt)
            if v_start_rel < 0:  # Clip left side of variant
                alt_length += v_start_rel
                v_start_rel = 0

            if v_end_rel > length:  # Clip right side of variant
                alt_length -= length - v_end_rel
                v_end_rel = length

            assert v_end_rel >= 0
            assert v_start_rel >= 0
            length += alt_length - (v_end_rel - v_start_rel)

        return length

    def __getstate__(self) -> bytes:
        return pickle.dumps([self.reference_genome, self.chromosome, self.start, self.ref, self.alt])

    def __setstate__(self, state: bytes) -> None:
        refg, chrom, start, ref, alt = pickle.loads(state)
        self.__init__(chrom, start, ref, alt, refg)


########################################################################


@_cxx.register
class VariantTable(_cxx.VariantTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all variants that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all variants that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all variants that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all variants that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all variants that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all variants that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that overlap `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all variants that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [v for v in variants if v.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
            All variants that span exactly `interval`.
        """
        mock_unreachable()
        return [Variant()]

    @mock
    @property
    def stranded(self):  # pragma: no cover
        """If `True` then the strand is significant when calling the `find_x` methods.

        Returns
        -------
        :py:class:`bool`
            Whether this table can contain negative stranded intervals.
        """
        return mock_result(bool)

    @mock
    def __getitem__(self, index):  # pragma: no cover
        """Access to all variants.

        Allows iteration over all variants in the table::

            for variant in table:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested junction.

        Returns
        -------
        :py:class:`~genome_kit.Variant`
           The variant at the given index.
        """
        return mock_result(Variant)

    @mock
    def where(self, mask):  # pragma: no cover
        """Filter variants by numpy mask.

        Fast extraction of table elements based on a mask::

            # Intended to be used like this.
            variants = table.where(mask)

            # It is faster but equivalent to this.
            variants = [table[i] for i in np.where(mask)[0]]

        Parameters
        ----------
        mask : :py:class:`ndarray`
            A boolean mask the same size as the table.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Variant`
           A list of variants selected by the mask.
        """
        mock_unreachable()
        return [Variant()]

    def __repr__(self):
        return "<VariantTable, len() = {}>".format(len(self))
