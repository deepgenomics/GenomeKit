# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import numpy as np

from . import _cxx
from ._cxx_util import mock, mock_result, mock_unreachable
from .interval import Interval

#########################################################################


@_cxx.register
class GenomeTrackBuilder(_cxx.GenomeTrackBuilder):
    """Builds a .gtrack file.

    Creating a .gtrack file with :py:class:`~genome_kit.GenomeTrackBuilder`
    generally involves four steps:

    1. Open the file for writing.
    2. Set further configuration options.
    3. Set data on intervals, either as numpy arrays or from disk.
    4. Call :py:meth:`~genome_kit.GenomeTrackBuilder.finalize`.

    For example, assuming we had a list of genomic intervals and
    corresponding data arrays we wanted to fill them with::

        # 1. Open a .gtrack 8-bit encoding of floats in range [0, 1]
        track = GenomeTrackBuilder("foo.gtrack", "f8", "single_stranded", Genome("hg19"))

        # 2. Configure the track further
        track.set_default_value(0)          # Fill gaps in the data with 0

        # 3. Set data for each interval
        for interval, data in entries:
            track.set_data(interval, data)

        # 4. Finish writing the file
        track.finalize()

    Tracks can have any dimension (the `dim` argument) or resolution
    (the `resolution` argument).
    Resolution allows track data to be specified and stored at coarser
    than 1bp resolution.

    For example, a track with `resolution=5` can
    load a WIG file with step/span of 5bp::

       fixedStep chrom=chr1 start=1 step=5 span=5
       0.57
       0.23
       ...

    The resulting track is equivalent to the following::

       0    1    2    3    4    5    6    7    8    9    ...  <-- position
       0.57 0.57 0.57 0.57 0.57 0.23 0.23 0.23 0.23 0.23 ...  <-- decoded value
    """

    __slots__ = ()  # <--- ADDING SLOTS IS OK

    @mock
    def __init__(self, outfile, etype, strandedness, reference_genome, dim=1, resolution=1):  # pragma: no cover
        """Open a .gtrack for writing.

        Parameters
        ----------
        outfile : :py:class:`str`
            The .gtrack file to create.
        etype : :py:class:`str`
            The encoding format.
            See :py:class:`~genome_kit.GenomeTrack` for details.
        strandedness : :py:class:`str`
            Determines the order by which the data is applied.

            * ``"single_stranded"``: both strands share the same data. The data is applied in Interval coordinate (reference strand) order.
            * ``"strand_unaware"``: ignores the Interval strand, data is applied in Interval coordinate (reference strand) order.
            * ``"strand_aware"``: data is applied from 5" end to 3" end (sense strand order).
        reference_genome : :class:`~genome_kit.Genome`
            The reference genome for the track data to build.
        dim : :py:class:`int`
            Optional. The number of dimensions in the resulting track.
        resolution : :py:class:`int`
            Optional. The resolution of the track data, in genomic positions.
        """

    @mock
    def set_default_value(self, value):  # pragma: no cover
        """Set the fill value for gaps in the track.

        The default value does not need to be encodable.
        For example, a default value of -1.0 for a dictionary that
        does not contain -1.0 is fine.

        Parameters
        ----------
        value : :py:class:`float`
            The value to fill gaps with.
        """
        mock_unreachable()

    @mock
    def set_sparsity(self, min_run=48, min_delta=0.0):  # pragma: no cover
        """Enable automatic sparsification of data intervals.

        Reduces file size of tracks that are dominated by `default_value`,
        with non-default values appearing at sparse intervals.
        *Highly recommended for variableStep WIG files.*

        The default behaviour of a track is to encode all data
        passed to :py:meth:`~genome_kit.GenomeTrackBuilder.set_data`.
        After calling :py:meth:`~genome_kit.GenomeTrackBuilder.set_sparsity`,
        the encoder will analyze the data arrays passed in and avoid
        encoding runs of `default_value`.

        Specifically, if there exists a run of at least `min_run`
        values such that each ``abs(data[i] - default_value) <= min_delta``,
        then that entire run will be treated as `default_value` and
        not explicitly stored in the file.

        Note that, in order to exclude a run from the track, a 12-byte entry
        must be addd to the track index, so setting `min_run` to very small values
        can actually increase file size and will definitely increase query time.

        Parameters
        ----------
        min_run : :py:class:`int`
            Optional. The minimum length of a run before it's excluded.
        min_delta : :py:class:`float`
            Optional. The minimum delta for a value to be considered distinct
            from `default_value`.
        """
        mock_unreachable()

    @mock
    def set_clamping(self):  # pragma: no cover
        """Clamp all values to the encodable range.

        By default, a value outside the encodable range will trigger
        an error, to notify the user that his/her data may not be encoded
        properly by the current format.

        Calling this method shows intent that values should be clamped
        to the encodable range. Currently allowed only for f-type encodings.
        """
        mock_unreachable()

    @mock
    def set_dict(self, dict):  # pragma: no cover
        """Set the dictionary for etypes `f2..8`.

        The dictionary must be a float16/32 numpy array.
        For etype f<n>, there must be `2**n` entries in the array.
        Numerical values must be in non-decreasing order.
        The last entry in `dict` can be :py:data:`np.nan`
        if you wish to be able to encode NaN at specific positions.

        Parameters
        ----------
        dict : :py:class:`np.ndarray`
            The array of values comprising the dictionary.
        """
        mock_unreachable()

    @mock
    def set_transform(self, a, b, c, d):  # pragma: no cover
        """Transform data from range `[a, b]` to range `[c, d]`.

        This specifies that any time data is about to be encoded, it
        will undergo an affine transformation from numeric range `[a, b]`
        to range `[c, d]`.

        Parameters
        ----------
        a,b : :py:class:`float`
            The source range.
        c,d : :py:class:`float`
            The target range.
        """
        mock_unreachable()

    @mock
    def set_restriction(self, restriction):  # pragma: no cover
        """Restrict data to a specific interval.

        The restriction interval acts as a allowed region, so that all
        data outside that region is either ignored or cropped away by
        :py:meth:`~genome_kit.GenomeTrackBuilder.set_data`.

        The restriction interval's strand is ignored, i.e. the restriction
        always allows data on either strand.

        This method is useful for making small version of full-sized data pipeline,
        for the sake of creating unit tests or iterating faster during development.

        If the restriction interval is not aligned with the track resolution, some
        then the restriction interval will be expanded until it is aligned.

        Parameters
        ----------
        restriction : :py:class:`~genome_kit.Interval`
            The interval on which to keep data.
        """
        mock_unreachable()

    @mock
    def set_data(self, interval, data):  # pragma: no cover
        """Set the track data on a specific interval.

        If the interval overlaps a previously-specified interval,
        an error is raised.

        The `dtype` of the data array must be compatible with the
        `etype` of the track being built.

        Normally, `data` must be an NxM array where N is
        the length of `interval` and M is `dim`.

        If track resolution R is greater than 1, then `data`
        must be (N/R)xM where N is the length of `interval`.
        The start and end of `interval` must be aligned to R.
        In other words, data is supplied (and stored) at
        a coarsened resolution.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The interval to fill with data.
        data : :py:class:`np.ndarray`
            The data array, ordered according to the strandedness argument passed to the
            builder's constructor.
        """
        mock_unreachable()

    @mock
    def set_data_from_wig(self, pos_strand_file, neg_strand_file=None):  # pragma: no cover
        """Load all data from a .wig or .wig.gz file.

        Loads all data from a fixedStep or variableStep WIG file.

        If the track is stranded, then the negative strand
        track data will be loaded from the second file.

        The number of data columns in the WIG must match the track dimension.
        The span of the data must match the track resolution, with the exception
        of the very last datum on a chromosome.

        Adjacent WIG intervals will automatically be merged into a single GTRACK
        interval, with different data encoded for each spanned position.
        Consider calling :py:meth:`~genome_kit.GenomeTrackBuilder.set_sparsity`.

        Parameters
        ----------
        pos_strand_file : :py:class:`str`
            Path to the .wig file for positive strand.
        neg_strand_file : :py:class:`str`
            Optional. Path to the .wig file for negative strand, if stranded.
        """
        mock_unreachable()

    @mock
    def set_data_from_bedgraph(self, pos_strand_file, neg_strand_file=None):  # pragma: no cover
        """Load all data from a .bedgraph or .bedgraph.gz file.

        Loads all data from a BEDGRAPH file.
        A direct analogue of :py:meth:`~genome_kit.GenomeTrackBuilder.set_data_from_wig`.

        Parameters
        ----------
        pos_strand_file : :py:class:`str`
            Path to the .bedgraph file for positive strand.
        neg_strand_file : :py:class:`str`
            Optional. Path to the .bedgraph file for negative strand, if stranded.
        """
        mock_unreachable()

    @mock
    def set_data_from_bed(self, bedfile, categories=None):  # pragma: no cover
        """Load all data from a .bed or .bed.gz file.

        Loads each interval from a BED file and incorporates it into the track.

        If `categories` is not specified, then the value contained in the `score`
        column of the BED file will be stored in the corresponding interval.

        If `categories` is specified, then the `name` column of the BED file
        determines the value stored on each interval. The interval's `name` will
        be used to find an index in `categories`, and that index will be stored.

        The strandedness of the track must match the BED file (i.e. an unstranded
        track expects the `strand` column to contain ``'.'``).
        The track dimension must be 1.
        The start and end of each interval must match the track resolution,
        with the exception of the very last position of a chromosome.

        Adjacent BED intervals will automatically be merged into a single contiguous
        GTRACK interval.
        In a GTRACK file, the data value of an interval is repeated for the size of the
        interval, unlike in a BED file where it's specified only once per interval.
        This can result in large GTRACK files, especially when at 1bp resolution (the default).
        Consider using :py:meth:`~genome_kit.GenomeTrackBuilder.set_default_value`
        and :py:meth:`~genome_kit.GenomeTrackBuilder.set_sparsity` to manage file size.

        Parameters
        ----------
        bedfile : :py:class:`str`
            Path to the .bed file.
        categories : :py:class:`list` of :py:class:`str`
            Optional. List of categories to match with the `name` column
            in the BED file.
        """
        mock_unreachable()

    @mock
    def flush(self):  # pragma: no cover
        """Flush the currently accumulated data to disk.

        This function can free up memory while building a track,
        but should only be called immediately after an entire chromosome
        has had its data set. Once this is called all chromosomes
        that currently have *any* data on them become fixed.
        """
        mock_unreachable()

    @mock
    def finalize(self):  # pragma: no cover
        """Finish writing the file and close it.

        This function must be called last, to complete creation of the GTRACK file.
        """
        mock_unreachable()

    @mock
    @property
    def dim(self):  # pragma: no cover
        """The dimensionality of the track.

        Returns
        -------
        :py:class:`int`
           The number of columns in each returned array.
        """
        return mock_result(int)

    @mock
    @property
    def resolution(self):  # pragma: no cover
        """The resolution of the track.

        Returns
        -------
        :py:class:`int`
           The number of genomic positions spanned by each datum.
        """
        return mock_result(int)

    @mock
    @property
    def stranded(self):  # pragma: no cover
        """The strandedness of the track.

        Returns
        -------
        :py:class:`bool`
           Whether the negative strand stores data unique from the positive strand.
        """
        return mock_result(bool)

    @mock
    @property
    def etype(self):  # pragma: no cover
        """The encoding type of the track.

        Returns
        -------
        :py:class:`str`
           A string identifying the encoding type.
        """
        return mock_result(str)

    @mock
    @property
    def dtype(self):  # pragma: no cover
        """The default dtype for the track.

        Returns
        -------
        :py:class:`type`
           The numpy dtype that this track will decode to by default.
        """
        return mock_result(type)

    @mock
    @property
    def reference_genome(self):  # pragma: no cover
        """The name of the reference genome the track will be mapped to.

        Returns
        -------
        :py:class:`str`
           The name of the reference genome, e.g. "hg38".
        """
        return mock_result(str)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to `.gtrack` file being written to.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/data/phastcons.gtrack"
        """
        return mock_result(str)

    @mock
    @property
    def data_size(self):  # pragma: no cover
        """The number of positions spanned by data in the GTRACK file.

        If sparsity is disabled, the data size is the number of positions
        spanned by data via
        :py:meth:`~genome_kit.GenomeTrackBuilder.set_data`.

        Enabling sparsity may decrease the data size.

        Using resolution > 1bp will also decrease the data size, since data is
        encoded into a coarse coordinate system, not at full genomic resolution.

        The number of bits that the data occupies in the GTRACK file is
        approximately `data_size * dim * bits_per_datum`, but the exact
        computation is non-trivial, so portion of GTRACK file dedicated
        to data may be higher or lower than this rough estimate.

        Returns
        -------
        :py:class:`int`
           The number of data spanned by encoded data.
        """
        return mock_result(int)

    @mock
    @property
    def index_size(self):  # pragma: no cover
        """The number of separate intervals in the GTRACK index.

        If sparsity is disabled, the index size is the number of times
        :py:meth:`~genome_kit.GenomeTrackBuilder.set_data` was called.

        Enabling sparsity may increase the index size.

        Returns
        -------
        :py:class:`int`
           The number of unique intervals in the index.
        """
        return mock_result(int)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __repr__(self):
        return "<GenomeTrackBuilder '{}'>".format(self.filename)


#########################################################################


@_cxx.register
class GenomeTrack(_cxx.GenomeTrack):
    """Access to a genomic data track.

    A :py:class:`~genome_kit.GenomeTrack` provides a read-only view of a `.gtrack` file.
    A track can be easily opened and queried::

        >>> interval = Interval('chr17', '+', 41246881, 41246886, 'hg19')
        >>> track = GenomeTrack("hg19.phastcons_mammal.f4.gtrack")
        >>> track(interval)
        array([[ 1.   ],
               [ 0.857],
               [ 0.571],
               [ 0.357],
               [ 0.286]], dtype=float16)

        >>> track(interval.as_opposite_strand())
        array([[ 0.286],
               [ 0.357],
               [ 0.571],
               [ 0.857],
               [ 1.   ]], dtype=float16)

    You can also open a track with a specific scope::

        >>> with GenomeTrack("foo.gtrack") as track:
        ...     values = track(interval)

    Tracks can be encoded in formats called `etypes`.
    The available etypes trade representation power
    for storage size and query speed::

        etypes = [  # DECODABLE TYPES          RANGE                ENCODED TYPE
            "m0",   # float16/32, uint8, bool  {1}                  0-bit (mask)
            "u1",   # float16/32, uint8, bool  {0..1}               1-bit
            "u2",   # float16/32, uint8        {0..3}               2-bit
            "u3",   # float16/32, uint8        {0..7}               3-bit
            "u4",   # float16/32, uint8        {0..15}              4-bit
            "u5",   # float16/32, uint8        {0..31}              5-bit
            "u6",   # float16/32, uint8        {0..63}              6-bit
            "u8",   # float16/32, uint8        {0..255}             8-bit
            "i8",   # float16/32, int8         {-128..127}          8-bit
            "f2",   # float16/32               [-65504,65504]       2-bit dictionary
            "f3",   # float16/32               [-65504,65504]       3-bit dictionary
            "f4",   # float16/32               [-65504,65504]       4-bit dictionary
            "f5",   # float16/32               [-65504,65504]       5-bit dictionary
            "f6",   # float16/32               [-65504,65504]       6-bit dictionary
            "f8",   # float16/32               [-65504,65504]       8-bit dictionary
            "f16",  # float16/32               [-65504,65504]       16-bit float
            "f32",  # float16/32               [-1.2e-38, 3.4e38]   32-bit float
        ]

    The default dictionary for `f2..8` types is the range [0,1] computed
    as ``np.linspace(0.0, 1.0, 2**num_bits, dtype=np.float16)``.

    Decoded types are called `dtypes`, and are indeed represented by
    numpy dtypes::

        dtypes = [
            np.bool_,    # {False, True}
            np.uint8,    # {0..255}
            np.int8,     # {-128..127}
            np.float16,  # [-65504, 65504] IEEE 754 half-precision
            np.float32,  # [-1.2e-38, 3.4e38]
        ]

    New tracks can be can created using :py:class:`~genome_kit.GenomeTrackBuilder`.
    """

    __slots__ = ()  # <--- ADDING SLOTS IS OK

    @mock
    def __init__(self, infile):  # pragma: no cover
        """Open a .gtrack file.

        Parameters
        ----------
        infile : :py:class:`str`
            The .gtrack file to open.
        """

    @mock
    def __call__(self, interval, dtype=None, out=None):  # pragma: no cover
        """Extract data from a genomic track.

        Data is always returned in sense-strand order, meaning it is
        always sensitive to the strand of the interval,
        even on an unstranded track.

        You can optionally specify the dtype that a track should be
        decoded to.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The stranded query interval.
        dtype : :py:class:`type`
            Optional. The numpy dtype of the resulting array.
        out : :py:class:`~numpy.ndarray`
            Optional. A numpy array to hold the result.

        Returns
        -------
        :py:class:`~numpy.ndarray`
            The decoded track data with ``len(interval)`` rows and ``dim`` columns.
        """
        return mock_result(np.ndarray)

    @mock
    @property
    def dim(self):  # pragma: no cover
        """The dimensionality of the track.

        Returns
        -------
        :py:class:`int`
           The number of columns in each returned array.
        """
        return mock_result(int)

    @mock
    @property
    def resolution(self):  # pragma: no cover
        """The resolution of the track.

        Returns
        -------
        :py:class:`int`
           The number of genomic positions spanned by each datum.
        """
        return mock_result(int)

    @mock
    @property
    def stranded(self):  # pragma: no cover
        """The strandedness of the track.

        Returns
        -------
        :py:class:`bool`
           Whether the negative strand stores data unique from the positive strand.
        """
        return mock_result(bool)

    @mock
    @property
    def etype(self):  # pragma: no cover
        """The encoding type of the track.

        Returns
        -------
        :py:class:`str`
           A string identifying the encoding type.
        """
        return mock_result(str)

    @mock
    @property
    def dtype(self):  # pragma: no cover
        """The default dtype for the track.

        Returns
        -------
        :py:class:`type`
           The numpy dtype that this track will decode to by default.
        """
        return mock_result(type)

    @mock
    @property
    def reference_genome(self):  # pragma: no cover
        """The name of the reference genome this track is mapped to.

        Returns
        -------
        :py:class:`str`
           The name of the reference genome, e.g. "hg38".
        """
        return mock_result(str)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to `.gtrack` file from which data is extracted.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/data/phastcons.gtrack"
        """
        return mock_result(str)

    @mock
    def intervals(self):  # pragma: no cover
        """The list of intervals covered by this track.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Interval`
           The list of intervals with data in this track.
        """
        return mock_result(str)

    def histogram(self, region, separate_dims=False):
        """Generate a histogram of the values in the track.

        Returns a dictionary where each key is a float appearing
        at least once in the track, and each value is a count of the
        number of times it occurred. NaN values are omitted.

        If `region` is a :py:class:`~genome_kit.Genome`, computes
        a genome-wide histogram.
        If `region` is an interval or list of intervals, computes the
        histogram over that region only.

        If `separate_dims` is true, returns a list of
        :py:attr:`~genome_kit.GenomeTrack.dim` histograms,
        one for each track dimension.

        Parameters
        ----------
        region : :py:class:`~genome_kit.Genome` | :py:class:`~genome_kit.Interval` | \
                :py:class:`list` of :py:class:`~genome_kit.Interval`
            Specific interval(s) upon which to compute
            the histogram.
        separate_dims : :py:class:`bool`
            Optional. Whether to generate a single histogram over all
            dimensions, or to return a separate histogram for each.

        Returns
        -------
        :py:class:`dict` | :py:class:`list` of :py:class:`dict`
            A dictionary ``{ value : count }`` where the `value`
            appearing in the track `count` times. Zero
        """
        # Get `intervals` list from the `interval` arg
        from .genome import Genome
        if isinstance(region, Genome):
            intervals = self._intervals_for_genome(region)
        elif isinstance(region, Interval):
            intervals = [region]
        elif isinstance(region, list):
            intervals = region
        else:
            raise TypeError("Expected region to be a Genome, Interval, or list of Interval")

        # Build an array where counters[i,j] is the number of times in dimension i that
        # the float16 value represented by j=0..2^16 appeared.
        counters = np.zeros((self.dim, 2**16), dtype=np.int64)
        for interval in intervals:
            counters += self._histogram_counters_for_interval(interval)

        # Collapse each dimension into one single histogram by summing up the counts unless separate_dims=True
        if separate_dims:
            return [self._histogram_dict_from_counters(row) for row in counters]
        else:
            return self._histogram_dict_from_counters(counters.sum(axis=0))

    def _intervals_for_genome(self, genome):
        """Returns a list of intervals, one for each chromosome, which effectively partition the entire reference genome
        """
        if genome.refg != self.refg:
            raise ValueError("Track reference genome '{}' does not match argument '{}'".format(self.refg, genome.refg))
        strands = ('+', '-') if self.stranded else ('+', )
        genome_intervals = []
        for chrom in genome.chromosomes:
            chrom_size = genome.chromosome_size(chrom)
            for strand in strands:
                genome_intervals.append(Interval(chrom, strand, 0, chrom_size, genome))

        return genome_intervals

    def _histogram_counters_for_interval(self, interval, slice_size=100000):
        """Returns an int64 array `counters` of size (dim, 2**16) where `counters[i,j]` is
        the number of times float16 value indexed by j appeared in dimension i in the
        given interval."""
        counters = np.zeros((self.dim, 2**16), dtype=np.int64)  # int64 for compatibility with np.bincount

        # Outer loop partitions each interval into chunks, to keep peak memory consistent
        for start in range(interval.start, interval.end, slice_size):
            end = min(start + slice_size, interval.end)

            # build an interval for this particular slice, and read the track values for it.
            interval_slice = Interval(interval.chromosome, interval.strand, start, end, interval.refg)
            slice_data_float16 = self(interval_slice, dtype=np.float16)
            try:
                slice_data_uint16 = np.frombuffer(
                    memoryview(slice_data_float16), dtype=np.uint16).reshape(slice_data_float16.shape)
            except AttributeError:
                slice_data_uint16 = np.frombuffer(
                    np.getbuffer(slice_data_float16), dtype=np.uint16).reshape(slice_data_float16.shape)

            # increment the appropriate key in counters (int_indexes[dim_index]) for each dimension
            # corresponding to the value of the track for that dimension
            for dim in range(self.dim):
                counters[dim] += np.bincount(slice_data_uint16[:, dim].astype(np.int64), minlength=2**16)

        return counters

    def _histogram_dict_from_counters(self, counters):
        """Converts a counters array into a dict, where each key is a float and each value is a non zero count.
        """
        # counters[i] is the number of times the float16 value float_keys[i]
        # appeared, for i=0..65536. Build a histogram of all non-nan keys.

        # merge the count for -0.0 into the count for +0.0
        pos_zero_index = 0  # index of +0.0
        neg_zero_index = 2**15  # index of -0.0
        counters[pos_zero_index] += counters[neg_zero_index]
        counters[neg_zero_index] = 0

        # generate key for all 65536 possible float16 values, so that
        # float_keys[i] is the float16 value counted by counters[i].
        int_keys = np.arange(2**16, dtype=np.uint16)
        try:
            float_keys = np.frombuffer(memoryview(int_keys), dtype=np.float16)
        except AttributeError:
            float_keys = np.frombuffer(np.getbuffer(int_keys), dtype=np.float16)

        # make dict of all non-nan keys and all non-zero values
        mask = np.logical_and(~np.isnan(float_keys), counters != 0)
        histogram = dict(zip(float_keys[mask], counters[mask]))

        return histogram

    @mock
    def close(self):  # pragma: no cover
        """Close the file handle.

        Note that if you use a `with` statement the file handle is automatically closed::

            # Open the file
            with GenomeTrack('foo.gtrack') as track:
                ...
            # <-- File is now closed
        """
        mock_unreachable()

    @mock
    @staticmethod
    def gtrack_version():  # pragma: no cover
        """The `gtrack` file format version that this build supports.

        Attempting to open a `gtrack` file of a mismatched version will raise an :py:exc:`IOError`.

        Returns
        -------
        :py:class:`int`
           The version number.
        """
        return mock_result(int)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __repr__(self):
        return "<GenomeTrack '{}'>".format(self.filename)
