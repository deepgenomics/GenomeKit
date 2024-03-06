# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from . import _cxx
from ._cxx_util import mock
from ._cxx_util import mock_result
from ._cxx_util import mock_unreachable
from ._cxx_util import strip_mock_bases
from .genome import Genome
from .interval import Interval
from .variant import Variant, VariantTable

#########################################################################


@_cxx.register
class ReadAlignments(_cxx.ReadAlignments):
    """Access to read alignments.
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def __init__(self, infile):  # pragma: no cover
        """Open a .ralign file.

        Parameters
        ----------
        infile : :py:class:`str`
            The .ralign file to open.

        Raises
        ------
        IOError
            Occurs when `infile` refers to an incompatible ralign file for this version of ``genome_kit``.
        """

    @mock
    @property
    def junctions(self):  # pragma: no cover
        """Access to junctions inferred from read alignments

        Returns
        -------
        :py:class:`~genome_kit.JunctionTable`
            A table with optimized methods to query over all junctions by position.
        """
        return mock_result(JunctionTable)

    @mock
    @property
    def alignments(self):  # pragma: no cover
        """Access to read alignments, i.e. a row in BAM file

        Returns
        -------
        :py:class:`~genome_kit.AlignmentTable`
            A table with optimized methods to query over all read alignments by position.
        """
        return mock_result(AlignmentTable)

    @mock
    @property
    def matches(self):  # pragma: no cover
        """Access to aligned regions of read alignments

        Returns
        -------
        :py:class:`~genome_kit.AlignmentMatchTable`
            A table with optimized methods to query over all aligned regions by position.
        """
        return mock_result(AlignmentMatchTable)

    @mock
    @property
    def variants(self):  # pragma: no cover
        """Access to variants observed by read alignments

        Returns
        -------
        :py:class:`~genome_kit.VariantTable`
            A table where each element (each row) is a variant that was observed
            in at least one read alignment.
        """
        return mock_result(VariantTable)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to the file from which read alignments are retrieved.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/mnt/data/sample000132.ralign"
        """
        return mock_result(str)

    @mock
    def close(self):  # pragma: no cover
        """Close the file handle.

        Note that if you use a `with` statement the file handle is automatically closed::

            # Open the file
            with ReadAlignments('my_file.ralign') as raligns:
                ...

            # <-- File is now closed
        """
        mock_unreachable()

    @mock
    @staticmethod
    def ralign_version():  # pragma: no cover
        """The `ralign` file format version that this build supports.

        Attempting to open an `ralign` file of a mismatched version will raise an :py:exc:`IOError`.

        Returns
        -------
        :py:class:`int`
           The version number.
        """
        return mock_result(int)

    @mock
    @staticmethod
    def build_ralign(
            outfile,
            infiles,
            reference_genome,
            exclude=None,
            allow=None,
            include_duplicates=False,
            library_format='U',
    ):  # pragma: no cover
        """Build an `ralign` file from one or more `SAM` files.

        If multiple input files are specified, their junctions and alignments will be pooled.

        Once an `ralign` file is created, it can be opened by creating
        a :py:class:`~genome_kit.ReadAlignments` object.

        **Alignments on the reference genome.**
        All junction, alignment, match, and variant intervals in the resulting file
        are on the reference genome.

        For example, an alignment starting at reference genome position 2
        and with CIGAR string ``1M1D1M2N1M``
        is matched to the reference genome as::

            0 1 2 3 4 5 6 7 8 9     <-- reference genome coordinates
                M D M - - M M       <-- M = match, D = deletion


        This results in one alignment (2:9), two matches (2:5, 7:9),
        one junction (5:7), and one variant (3:4).

        Similarly for insertions, with CIGAR string ``1M1I1M2N1M`` aligns as::

            0 1 2 3 4 5 6 7 8 9
                M M - - M M
                 ^
                 I                  <-- I = insertion between

        This results in one alignment (2:8), two matches (2:4, 6:8),
        one junction (4:6), and one variant (3:3).


        Parameters
        ----------
        outfile : :py:class:`str`
            The path to the destination `rdist` file.

        infiles : :py:class:`list` of :py:class:`str` | :py:class:`file`
            The paths to the source `SAM` files, or a reference to ``sys.stdin``.
            Streaming lines from stdin is useful for reading directly from BAM
            via ``samtools view -h``. Be sure to use ``-h`` so that
            the reference genome can be inferred from the header lines.

        reference_genome : :class:`~genome_kit.Genome`
            The reference genome of the data in ``infiles``.
            If the genome contains annotations, junction strands will be
            inferred and validated against the introns annotations.

        exclude : :py:class:`list` of :py:class:`~genome_kit.Interval`, optional
            Junctions within these intervals will be excluded.

        allow : :py:class:`list` of :py:class:`~genome_kit.Interval`, optional
            Only junctions within these intervals will be included,
            so long as they are not excluded.

        include_duplicates : :py:class:`bool`, optional
            Includes duplicates when `True`. Defaults to `False`. See
            https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf
            subsection 1.4.2

        library_format : :class:`str`, optional
            Specifies to infer junction strands via the library format string.
            "SF" implies the first read is the sense strand, whereas "SR"
            implies the second read is the sense strand. (the default is 'U',
            which represents unstranded [no strand inference]).
            See https://salmon.readthedocs.io/en/stable/library_type.html.
        """
        mock_unreachable()

    PILEUP_A = 0
    PILEUP_C = 1
    PILEUP_G = 2
    PILEUP_T = 3
    PILEUP_DEL = 4

    def pileup(self, interval, dtype=np.int32):  # pragma: no cover
        """
        Returns the DNA sequence pileup for body reads in this `ReadAlignments`
        as a 5-track. The pileup is variant-aware, so the track counts
        the alternate sequences, instead of the reference, on any alignment
        matches with variants.

        Currently, only acgtACGT substitutions and deletions variants are
        supported.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            A positive interval to query over all the read alignments.
        dtype : :py:class:`~numpy.dtype`, optional
            The data type for the returned counts.

        Returns
        -------
        :py:class:`~numpy.ndarray`
            An array of ``shape`` (len(`interval`), 5) where the second
            dimension is indexed via ``PILEUP_A``, ``PILEUP_C``, ``PILEUP_G``,
            ``PILEUP_T``, and ``PILEUP_DEL``.

        """
        # "matches" refers to aligned regions: this includes mismatches.
        # increment all the aligned regions and then fix-up for variant regions
        counts = np.zeros((len(interval), 5), dtype)
        ref_counts = np.zeros(len(interval), dtype)

        for match in self.matches.find_overlapping(interval):
            start = max(match.start, interval.start) - interval.start
            end = min(match.end, interval.end) - interval.start
            ref_counts[start:end] += 1

            for variant in (x for x in match.variants if x.overlaps(interval.as_positive_strand())):
                # ignore inserts for now
                diff_len = len(variant.alt) - len(variant.ref)
                if diff_len > 0:
                    continue

                start = max(variant.start, interval.start)
                end = min(variant.end, interval.end)
                count_start = start - interval.start
                count_end = end - interval.start
                ref_counts[count_start:count_end] -= 1
                position_indices = np.arange(count_start, count_end)

                if diff_len == 0:
                    alt = variant.alt[start - variant.start:end - variant.start]
                    n_indices = [i for i, x in enumerate(alt) if x == 'N']
                    pileup_indices = np.delete(self._dna_as_pileup_index(alt), n_indices)
                    counts[np.delete(position_indices, n_indices), pileup_indices] += 1
                elif diff_len < 0:
                    assert not variant.alt  # normalized to only include deleted region (no padding)
                    counts[position_indices, self.PILEUP_DEL] += 1

        assert (not ref_counts[ref_counts < 0].any())
        dna_str = Genome(interval.reference_genome).dna(interval.as_positive_strand())
        n_indices = [i for i, x in enumerate(dna_str) if x == 'N']
        position_indices = np.delete(np.arange(len(interval)), n_indices)
        pileup_indices = np.delete(self._dna_as_pileup_index(dna_str), n_indices)
        ref_counts = np.delete(ref_counts, n_indices)
        counts[position_indices, pileup_indices] += ref_counts
        return counts

    DNA_TO_INDEX = [b'\xFF'] * 256
    DNA_TO_INDEX[ord('a')] = DNA_TO_INDEX[ord('A')] = b'\0'
    DNA_TO_INDEX[ord('c')] = DNA_TO_INDEX[ord('C')] = b'\1'
    DNA_TO_INDEX[ord('g')] = DNA_TO_INDEX[ord('G')] = b'\2'
    DNA_TO_INDEX[ord('t')] = DNA_TO_INDEX[ord('T')] = b'\3'
    DNA_TO_INDEX = b"".join(DNA_TO_INDEX)

    @classmethod
    def _dna_as_pileup_index(cls, dna_str):
        """
        Convert a case-insensitive DNA string (ACGT) into an array of uint8 indices
        in range {0,1,2,3}. Non-ACGT characters will be converted to index 255.
        """
        indices = dna_str.encode('ascii').translate(cls.DNA_TO_INDEX)
        return np.ndarray(shape=(len(indices), ), buffer=indices, dtype=np.uint8)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __repr__(self):
        return '<ReadAlignments "{}">'.format(self.filename)


########################################################################
#   Junction, JunctionTable
########################################################################


@strip_mock_bases
@_cxx.register
class Junction(_cxx.Junction, Interval):
    """A Junction represents a gap in a read alignment when compared to the reference genome.

    Bases: :py:class:`~genome_kit.Interval`
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    # noinspection PyMissingConstructor
    @mock
    def __init__(self):  # pragma: no cover
        pass  # Stub to prevent Interval.__init__ from being accessible to docs

    @mock
    @property
    def interval(self):  # pragma: no cover
        """The interval spanned by this junction.

        Note that :py:class:`~genome_kit.Junction` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this junction.
        """
        return mock_result(Interval)

    @mock
    @property
    def _num_alignments(self):  # pragma: no cover
        """
        Returns
        -------
        :py:class:`int`
        """
        return mock_result(int)

    @mock
    def _get_alignment(self, index):  # pragma: no cover
        """
        Parameters
        ----------
        index : :py:class:`int`

        Returns
        -------
        :py:class:`~genome_kit.Alignment`
        """
        return mock_result(Alignment)

    class _AlignmentsTupleAccessor(object):
        __slots__ = '_junction'

        def __init__(self, junction):  # pragma: no cover
            self._junction = junction

        def __len__(self):  # pragma: no cover
            return self._junction._num_alignments

        def __getitem__(self, item):  # pragma: no cover
            return self._junction._get_alignment(item)

    @property
    def alignments(self):  # pragma: no cover
        """The read alignments across this junction.

        Returns
        -------
        :py:class:`tuple` of :py:class:`~genome_kit.Alignment`
            The read alignments across this junction.
        """
        return self._AlignmentsTupleAccessor(self)

    def __repr__(self):
        return '<Junction {}, len(alignments) = {}>'.format(self.interval, self._num_alignments)


@_cxx.register
class JunctionTable(_cxx.JunctionTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all junctions that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end5 == interval.end5]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junctions that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Junction()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all junctions that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end3 == interval.end3]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junctions that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Junction()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all junctions that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end5.expand(0, 1) in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junctions that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Junction()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all junctions that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end3.expand(1, 0) in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junctions that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Junction()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all junctions that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junction that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Junction()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all junctions that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.overlaps(interval)]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junctions that overlap `interval`.
        """
        mock_unreachable()
        return [Junction()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all junctions that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.interval == interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Junction`
            All junctions that span exactly `interval`.
        """
        mock_unreachable()
        return [Junction()]

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
        """Lookup a junction by index

        Allows iteration over all junctions::

            for junc in junctions:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested junction.

        Returns
        -------
        :py:class:`~genome_kit.Junction`
           The junction object identified by the given index.
        """
        return mock_result(Junction)

    def __repr__(self):
        return "<JunctionTable, len = {}>".format(len(self))


########################################################################
#   Alignment, AlignmentTable
########################################################################


@strip_mock_bases
@_cxx.register
class Alignment(_cxx.Alignment, Interval):
    """A read alignment

    Bases: :py:class:`~genome_kit.Interval`
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    # noinspection PyMissingConstructor
    @mock
    def __init__(self):  # pragma: no cover
        pass  # Stub to prevent Interval.__init__ from being accessible to docs

    @mock
    @property
    def interval(self):  # pragma: no cover
        """The interval spanned by this read alignment.

        Note that :py:class:`~genome_kit.Alignment` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this read alignment.
        """
        return mock_result(Interval)

    @mock
    @property
    def matches(self):  # pragma: no cover
        """The matching regions of this read alignment

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            A list of matching regions
        """
        mock_unreachable()
        return [AlignmentMatch()]

    def __repr__(self):
        return '<Alignment {}, len(matches) = {}>'.format(self.interval, len(self.matches))


@_cxx.register
class AlignmentTable(_cxx.AlignmentTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all read alignments that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read.end5 == interval.end5]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignments that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Alignment()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all read alignments that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read.end3 == interval.end3]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignments that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Alignment()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all read alignments that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read.end5.expand(0, 1) in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignments that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Alignment()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all read alignments that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read.end3.expand(1, 0) in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignments that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Alignment()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all read alignments that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignment that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Alignment()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all read alignments that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read.overlaps(interval)]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignments that overlap `interval`.
        """
        mock_unreachable()
        return [Alignment()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all read alignments that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [read for read in reads if read.interval == interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Alignment`
            All read alignments that span exactly `interval`.
        """
        mock_unreachable()
        return [Alignment()]

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
        """Lookup a read alignment by index

        Allows iteration over all read alignments::

            for read in reads:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested read alignment.

        Returns
        -------
        :py:class:`~genome_kit.Alignment`
           The read alignment object identified by the given index.
        """
        return mock_result(Alignment)

    def __repr__(self):
        return "<AlignmentTable, len = {}>".format(len(self))


########################################################################
#   AlignmentMatch, AlignmentMatchTable
########################################################################


@strip_mock_bases
@_cxx.register
class AlignmentMatch(_cxx.AlignmentMatch, Interval):
    """A matching region a read alignment

    Bases: :py:class:`~genome_kit.Interval`
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    # noinspection PyMissingConstructor
    @mock
    def __init__(self):  # pragma: no cover
        pass  # Stub to prevent Interval.__init__ from being accessible to docs

    @mock
    @property
    def interval(self):  # pragma: no cover
        """The interval spanned by the matching region of a read alignment.

        Note that :py:class:`~genome_kit.AlignmentMatch` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this alignment match.
        """
        return mock_result(Interval)

    @mock
    @property
    def variants(self):  # pragma: no cover
        """A list of variants observed in this matching region of an alignment.

        Returns
        -------
        :py:class:`tuple` of :py:class:`~genome_kit.Variant`
            A tuple of variants that this matching region of an alignment contains
        """
        mock_unreachable()
        return (Variant(), )

    def __repr__(self):
        return '<AlignmentMatch {}, len(variants) = {}>'.format(self.interval, len(self.variants))


@_cxx.register
class AlignmentMatchTable(_cxx.AlignmentMatchTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all alignment matches that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match.end5 == interval.end5]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignment matches that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all alignment matches that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match.end3 == interval.end3]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignment matches that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all alignment matches that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match.end5.expand(0, 1) in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignment matches that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all alignment matches that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match.end3.expand(1, 0) in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignments matches that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all alignment matches that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match in interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignment matches that fall entirely within `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all alignment matches that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match.overlaps(interval)]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignment matches that overlap `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all alignment matches that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [match for match in matches if match.interval == interval]

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.AlignmentMatch`
            All alignment matches that span exactly `interval`.
        """
        mock_unreachable()
        return [AlignmentMatch()]

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
        """Lookup an alignment match by index

        Allows iteration over all alignment matches::

            for match in matches:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested alignment match.

        Returns
        -------
        :py:class:`~genome_kit.AlignmentMatch`
           The alignment match object identified by the given index.
        """
        return mock_result(AlignmentMatch)

    def __repr__(self):
        return "<AlignmentMatchTable, len = {}>".format(len(self))
