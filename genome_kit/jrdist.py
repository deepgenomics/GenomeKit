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
import os

#########################################################################


@_cxx.register
class JunctionReadCount(_cxx.JRCount):
    """A junction read count at a specific shift from the junction.
    """
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    @property
    def strand(self):  # pragma: no cover
        """The strand that this count refers to.
        Returns
        -------
        :py:class:`str`
            The strand that this count refers to.
        """
        return mock_result(str)

    @mock
    @property
    def shift(self):  # pragma: no cover
        """The number of aligned positions between the far end (start/ end) of the alignment and the junction.

        If negative, the count refers to an alignment upstream of the junction interval
        (always defined on the positive strand). Does not count insertions, deletions, clippings, etc.

        Returns
        -------
        :py:class:`int`
            The shift relative to the junction.
        """
        return mock_result(int)

    @mock
    @property
    def count(self):  # pragma: no cover
        """The number of reads at the given strand and shift from junction.
        Returns
        -------
        :py:class:`int`
            The number of reads merged into the count.
        """
        return mock_result(int)

    def __repr__(self):
        return '<JunctionReadCount ({}, {}, {})>'.format(self.strand, self.shift, self.count)


########################################################################


@strip_mock_bases
@_cxx.register
class JunctionReadDistribution(_cxx.JRDist, Interval):
    """A read distribution across a junction.
    A junction read distribution comprises a list of non-zero counts across the junction,
    where the shift and strand of each count is preserved.

    Each read mapping across the junction corresponds to two counts:
    one for the part mapping upstream of the junction and another for the part mapping downstream of it.

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

        Note that :py:class:`~genome_kit.JunctionReadDistribution` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this junction.
        """
        return mock_result(Interval)

    @mock
    @property
    def num_counts(self):  # pragma: no cover
        """The number of counts on this junction.

        Returns
        -------
        :py:class:`int`
            The number of unique non-zero counts stored on this junction, equivalent to the `len` of this object.
        """
        return mock_result(int)

    @mock
    @property
    def num_reads(self):  # pragma: no cover
        """The number of reads across this junction.

        Returns
        -------
        :py:class:`int`
            The number of reads across this junction, equal to the summing all counts and dividing by two.
        """
        return mock_result(int)

    @mock
    def __getitem__(self, index):  # pragma: no cover
        """Access to a junction read count.

        Allows iteration over all read counts across this junction::

            for count in junction:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested count entry.

        Returns
        -------
        :py:class:`~genome_kit.JunctionReadCount`
           The junction read count identified by the given index.
        """
        return mock_result(JunctionReadCount)

    def __repr__(self):
        return '<JunctionReadDistribution {} ({} counts)>'.format(self.as_ucsc(), self.num_counts)


########################################################################


@_cxx.register
class JunctionReadDistributionTable(_cxx.JRDistTable):
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

        The results will also be sorted by 5'.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that fall entirely within `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that overlap `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`list` of :py:class:`~genome_kit.JunctionReadDistribution`
            All junctions that span exactly `interval`.
        """
        mock_unreachable()
        return [JunctionReadDistribution()]

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
        :py:class:`~genome_kit.JunctionReadDistribution`
           The junction read distribution object identified by the given index.
        """
        return mock_result(JunctionReadDistribution)

    def __repr__(self):
        return "<JunctionReadDistributionTable, len() = {}>".format(len(self))


########################################################################


@_cxx.register
class ReadDistributions(_cxx.ReadDistributions):
    """Access to read distributions.
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def __init__(self, infile):  # pragma: no cover
        """Open a .rdist file.

        Parameters
        ----------
        infile : :py:class:`str`
            The .rdist file to open.
        """

    @mock
    @property
    def junctions(self):  # pragma: no cover
        """Access to read distributions that are split across junctions.

        Returns
        -------
        :py:class:`~genome_kit.JunctionReadDistributionTable`
           An object that supports looping or queries over junctions, where
           each junction provides a read distribution.
        """
        return mock_result(JunctionReadDistributionTable)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to the file from which read distributions are retrieved.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/mnt/data/sample000132.rdist"
        """
        return mock_result(str)

    @mock
    def close(self):  # pragma: no cover
        """Close the file handle.

        Note that if you use a `with` statement the file handle is automatically closed::

            # Open the file
            with ReadDistributions('my_file.ralign') as raligns:
                ...

            # <-- File is now closed
        """
        mock_unreachable()

    @mock
    @staticmethod
    def rdist_version():  # pragma: no cover
        """The `rdist` file format version that this build supports.

        Attempting to open an `rdist` file of a mismatched version will raise an :py:exc:`IOError`.

        Returns
        -------
        :py:class:`int`
           The version number.
        """
        return mock_result(int)

    @mock
    @staticmethod
    def build_rdist(outfile,
                    infiles,
                    min_reads=0,
                    min_overhang=0,
                    exclude=None,
                    allow=None,
                    pass_filter=None):  # pragma: no cover
        """Build an `rdist` file from one or more `jralign` files.

        If multiple input files are specified, their read counts will be pooled.

        Once an `rdist` file is created, it can be opened by creating
        a :py:class:`~genome_kit.ReadDistributions` object .

        Parameters
        ----------
        outfile : :py:class:`str`
            The path to the destination `rdist` file.

        infiles : :py:class:`list` of :py:class:`str`
            The paths to the source `jralign` files.

        min_reads : :py:class:`int`
            Optional. The minimum number of reads across a junction to
            be considered for inclusion in the output file.

        min_overhang : :py:class:`int`
            Optional. The minimum overhang for a read to be included in a count.

        exclude : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Junctions within these intervals will be excluded.

        allow : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Only junctions within these intervals will be included,
            so long as they are not excluded.

        pass_filter : Callable[[str, int, int], bool], optional
            A pass filter to be applied while iterating over the
            (infile, junction_index, read_index). For example::

                # downsample to 10%
                pass_filter=lambda _, _, _: np.random.binomial(1, 0.1) != 0

                # only take reads that have variants
                jraligns = {x: JReadAlignments(x) for x in infiles}
                pass_filter=lambda infile, junction, read: jraligns[infile].junctions[junction][read].num_variants != 0
                ...
                for x in jraligns.values():
                    x.close()

        """
        mock_unreachable()

    @staticmethod
    def from_jraligns(infiles,
                      min_reads=0,
                      min_overhang=0,
                      exclude=[],
                      allow=[],
                      pass_filter=None,
                      cache=True,
                      outfile=None):
        """
        Returns a :class:`~genome_kit.ReadDistributions` object opened from a
        set of jralign files.

        This convenience function builds and then returns ReadDistributions in
        one step, by default the binary file (`.rdist`) is adjacent to the
        first input file.

        It is intended to be used as follows::

            rdist = ReadDistributions.from_jraligns(["input.jralign"])
            junctions = rdist.find_within(brca1_interval)
            ...

        If any of the filtering parameters or `infiles` contents change,
        `cache` should be set to `False` to regenerate the `ReadDistribution`.

        See Also
        --------
        :meth:`~genome_kit.ReadDistributions.build_rdist`

        Parameters
        ----------
        infiles : :py:class:`list` of :py:class:`str`
            The paths to the source `jralign` files.

        min_reads : :py:class:`int`
            Optional. The minimum number of reads across a junction to
            be considered for inclusion in the output file.

        min_overhang : :py:class:`int`
            Optional. The minimum overhang for a read to be included in a count.

        exclude : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Junctions within these intervals will be excluded.

        allow : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Only junctions within these intervals will be included,
            so long as they are not excluded.

        pass_filter : Callable[[str, int, int], bool], optional
            A pass filter to be applied while iterating over the
            (infile, junction_index, read_index). For example::

                # downsample to 10%
                pass_filter=lambda _, _, _: np.random.binomial(1, 0.1) != 0

                # only take reads that have variants
                jraligns = {x: JReadAlignments(x) for x in infiles}
                pass_filter=lambda infile, junction, read: jraligns[infile].junctions[junction][read].num_variants != 0
                ...
                for x in jraligns.values():
                    x.close()

        cache : :py:class:`bool`, optional
            Whether to re-use the binary or to re-generate it from `infiles`.

        outfile : :py:class:`str`, optional
            The path to the destination `rdist` file.

        Returns
        -------
        :class:`~genome_kit.ReadDistributions`
            The opened `rdist` file.
        """
        # If user doesn't specify path for binary output file, write it adjacent to the first jralign
        if outfile is None:
            outfile = os.path.splitext(infiles[0])[0] + ".rdist"  # pragma: no cover

        if cache and os.path.exists(outfile):
            if any(os.path.getmtime(x) > os.path.getmtime(outfile) for x in infiles):
                print("Warning: %s older than %s in ReadDistributions.from_jralign." %
                      (outfile, infiles))  # pragma: no cover
        else:
            ReadDistributions.build_rdist(outfile,
                                          infiles,
                                          min_reads=min_reads,
                                          min_overhang=min_overhang,
                                          exclude=exclude,
                                          allow=allow,
                                          pass_filter=pass_filter)

        # Open the binary file
        return ReadDistributions(outfile)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __repr__(self):
        return '<ReadDistributions "{}">'.format(self.filename)


########################################################################
