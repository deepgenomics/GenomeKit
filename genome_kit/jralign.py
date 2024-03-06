# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from . import _cxx
from ._cxx_util import mock
from ._cxx_util import mock_result
from ._cxx_util import mock_unreachable
from ._cxx_util import strip_mock_bases
from .interval import Interval
from .variant import Variant, VariantTable

#########################################################################


@_cxx.register
class JunctionReadAlignment(_cxx.JRAlign):
    """A junction read alignment.

    The `left` and `right` overhangs are positive integers representing the number
    of matched positions on either side of the junction.

    Overhangs are counted on the private genome from which the data was sequenced,
    and may not map to reference genome coordinates.
    As such, one should not expect meaningful reference genome coordinates by
    adding overhangs to the junction start/end.

    If `strand` is positive, then the read started `left` positions before the
    first intronic position on the private genome.
    If `strand` is negative, then the read started
    `right` positions after the last intronic position on the private genome.
    """
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    @property
    def left(self):  # pragma: no cover
        """The left overhang relative to the junction.

        Returns
        -------
        :py:class:`int`
            The number of matched positions left of the junction.
        """
        return mock_result(int)

    @mock
    @property
    def right(self):  # pragma: no cover
        """The right overhang relative to the junction.

        Returns
        -------
        :py:class:`int`
            The number of matched positions right of the junction.
        """
        return mock_result(int)

    @mock
    @property
    def strand(self):  # pragma: no cover
        """The strand ("+" or "-") that this read was aligned to.

        Not to be confused with the strand of the junction itself, which may or
        may not have been inferred.

        Returns
        -------
        :py:class:`chr`
            The strand the read was aligned to by the original aligner.
        """
        return mock_result(chr)

    @mock
    @property
    def num_variants(self):  # pragma: no cover
        """The number of variants observed in this read alignment.

        Returns
        -------
        :py:class:`int`
            The number of variants in the alignment.
        """
        return mock_result(int)

    @mock
    def variants(self):  # pragma: no cover
        """A list of variants observed in this read alignment.

        Returns
        -------
        :py:class:`tuple` of :py:class:`~genome_kit.Variant`
            A tuple of variants that this read contains
        """
        mock_unreachable()
        return (Variant(), )

    def __repr__(self):
        return '<JunctionReadAlignment ({}, {}, {}, {})>'.format(self.strand, self.left, self.right, self.num_variants)


########################################################################


@strip_mock_bases
@_cxx.register
class JunctionReadAlignments(_cxx.JRAligns, Interval):
    """A set of read alignments across a junction.

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

        Note: All junction intervals are on the reference strand, even if they comprise reads that were mapped to the
        reverse strand.
        When querying junctions, your query intervals must therefore (currently) be on the positive strand.

        Note that :py:class:`~genome_kit.JunctionReadAlignments` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this junction.
        """
        return mock_result(Interval)

    @mock
    @property
    def num_reads(self):  # pragma: no cover
        """The number of reads across this junction.

        Returns
        -------
        :py:class:`int`
            The number of reads across this junction, equivalent to the `len` of this object.
        """
        return mock_result(int)

    @mock
    def __getitem__(self, index):  # pragma: no cover
        """Access to a junction read alignment.

        Allows iteration over all read alignments across this junction::

            for align in junction:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested alignment entry.

        Returns
        -------
        :py:class:`~genome_kit.JunctionReadAlignment`
           The junction read alignment identified by the given index.
        """
        return mock_result(JunctionReadAlignment)

    def __repr__(self):
        return '<JunctionReadAlignments {} ({} reads)>'.format(self.as_ucsc(), self.num_reads)


########################################################################


@_cxx.register
class JunctionReadAlignmentsTable(_cxx.JRAlignsTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all junctions that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all junctions that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all junctions that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all junctions that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all junctions that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that fall entirely within `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all junctions that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that overlap `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all junctions that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [junc for junc in junctions if junc.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadAlignments`
            All junctions that span exactly `interval`.
        """
        mock_unreachable()
        return [JunctionReadAlignments()]

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
        """Access to all junctions.

        Allows iteration over all junctions in the read distribution table::

            for junc in junctions:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested junction.

        Returns
        -------
        :py:class:`~genome_kit.JunctionReadAlignments`
           The junction read alignment object identified by the given index.
        """
        return mock_result(JunctionReadAlignments)

    def __repr__(self):
        return "<JunctionReadAlignmentsTable, len() = {}>".format(len(self))


########################################################################


@_cxx.register
class JReadAlignments(_cxx.JReadAlignments):
    """Access to read alignments.
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def __init__(self, infile):  # pragma: no cover
        """Open a .jralign file.

        Parameters
        ----------
        infile : :py:class:`str`
            The .jralign file to open.
        """

    @mock
    @property
    def junctions(self):  # pragma: no cover
        """Access to read alignments that are split across junctions.

        Returns
        -------
        :py:class:`~genome_kit.JunctionReadAlignmentsTable`
             An object that supports looping or queries over junctions, where
             each junction provides a read distribution.

        Notes
        -----
        The returned table. junctions, and child objects should not be accessed
        after this JReadAlignments has been closed via :meth:`~.close`.
        """
        return mock_result(JunctionReadAlignmentsTable)

    @mock
    @property
    def variants(self):  # pragma: no cover
        """Access to variants observed over all read alignments

        Returns
        -------
        :py:class:`~genome_kit.VariantTable`
            A table where each element (each row) is a variant that was observed
            in at least one read alignment.

        Notes
        -----
        The returned table and variants should not be accessed after this
        JReadAlignments has been closed via :meth:`~.close`.
        """
        return mock_result(VariantTable)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to the file from which read alignments are retrieved.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/mnt/data/sample000132.jralign"
        """
        return mock_result(str)

    @mock
    def close(self):  # pragma: no cover
        """Close the file handle.

        Note that if you use a `with` statement the file handle is automatically closed::

            # Open the file
            with JReadAlignments('my_file.jralign') as jraligns:
                ...

            # <-- File is now closed

        Notes
        -----
        Closing this JReadAlignments will invalidate all tables and objects previously accessed via
        :meth:`~.junctions` and :meth:`~.variants`.
        """
        mock_unreachable()

    @mock
    @staticmethod
    def jralign_version():  # pragma: no cover
        """The `jralign` file format version that this build supports.

        Attempting to open an `jralign` file of a mismatched version will raise an :py:exc:`IOError`.

        Returns
        -------
        :py:class:`int`
           The version number.
        """
        return mock_result(int)

    @mock
    @staticmethod
    def build_jralign(outfile,
                      infiles,
                      reference_genome,
                      min_reads=0,
                      min_overhang=0,
                      exclude=None,
                      allow=None,
                      include_variants=False,
                      include_duplicates=False,
                      library_format='U',
                      overhang_error='error'):  # pragma: no cover
        """Build an `jralign` file from one or more `SAM` files.

        If multiple input files are specified, their junctions and read alignments will be pooled.

        Once an `jralign` file is created, it can be opened by creating
        a :py:class:`~genome_kit.JReadAlignments` object.

        Below are some details on how SAM alignments are processed, which are
        important for understanding how a CIGAR strings containing insertions/deletions
        (`5M8D5M`) or multiple splits (`5M100N10M50N5M`) are processed:

        * **Multi-split alignments.**
          If a read maps across multiple junctions, its CIGAR string will
          contain multiple ``N`` blocks. The read will then contribute a
          separate (left, right) alignment to *each* junction that it spans.

          For example, consider ``1S2M2N1M1N1M`` where the first match is
          aligned to start position 2::

              0 1 2 3 4 5 6 7 8 9
                S M M - - M - M

          This read will contribute a separate alignment to each junction.
          The `left` and `right` overhang is the sum of all alignment matches
          before/after that particular junction.
          In the above example, the resulting alignments would be
          be as if two separate reads were aligned as follows::

              0 1 2 3 4 5 6 7 8 9
                S M M - - M M        <-- alignment for 1st junction
                    S M M M - M      <-- alignment for 2nd junction

          In this scenario, the read would not actually match the reference genome
          in the above alignments, but this scheme facilitates counting of relative
          junction usage (assuming usage of the two junctions is independent)
          and also is allows positional bootstrap to be applied correctly
          without modification.

        * **Overhang filter.**
          The ``min_overhang`` filter is applied on a per-junction bases.
          For example, if ``min_overhang=5`` then a read mapped across
          two junctions ``20M100N10M100N3M`` will contribute an alignment
          only to the first, because the second has a right overhang of 3.

        * **Counting insertions/deletions.**
          Only alignment matches are counted. Other CIGAR types (e.g.,
          insertions, deletions, clippings) are not relevant when filtering
          for alignment quality.

          For example, a read starting at reference genome position 2
          and with CIGAR string ``1M1D1M2N1M``
          is matched to the reference genome as::

              0 1 2 3 4 5 6 7 8 9     <-- reference genome coordinates
                  M D M - - M         <-- M = match, D = deletion

          so the start of the junction is 2+1+1+1 = 5, but the
          the left overhang is only 2, because the deleted base
          doesn't exist in the private genome.

          Similarly for insertions, with CIGAR string ``1M1I1M2N1M``::

              0 1 2 3 4 5 6 7 8 9
                  M M - - M
                   ^
                   I                <-- I = insertion between

           the start of the junction is 2+1+1 = 4, but the left
           overhang is actually 2, because the inserted base exists
           only in the private genome.

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

        min_reads : :py:class:`int`
            Optional. The minimum number of reads across a junction to
            be considered for inclusion in the output file.

        min_overhang : :py:class:`int`
            Optional. The minimum overhang for a read to be included.

        exclude : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Junctions within these intervals will be excluded.

        allow : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Only junctions within these intervals will be included,
            so long as they are not excluded.

        include_variants : :py:class`bool`
            Optional. True if variants are to be included in the output file.
            Default will exclude all variants.

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

        overhang_error: : :class:`str`
            Determines how build_jralign handles reads with >255 length
            overhangs. "error" stops and reports the error. "clamp" limits
            the read to the maximum 255 positions on each side.
            This only affects overhangs and will not clip or filter variants.
        """
        mock_unreachable()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __repr__(self):
        return '<JReadAlignments "{}">'.format(self.filename)


########################################################################
