# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from . import _cxx
from ._cxx_util import mock
from ._cxx_util import mock_result

#########################################################################


@_cxx.register
class GenomeDNA(_cxx.GenomeDNA):
    """Access to the DNA sequence of a reference genome.

    This object should be accessed through the `dna` property of
    :py:class:`~genome_kit.Genome`::

        >>> from genome_kit import Genome, Interval
        >>> interval = Interval("chr5", "+", 50000, 50010, "hg19")
        >>> genome = Genome("hg19")
        >>> genome.dna(interval)
        'TAAACCACAT'

    """

    __slots__ = ()  # <--- ADDING SLOTS IS OK

    @mock
    def __call__(self, interval, allow_outside_chromosome=True):  # pragma: no cover
        """Extract DNA from a reference genome.

        If `interval` is on the negative strand, the reverse complement sequence is returned.
        If `interval` is outside the range of the chromosome, an :py:exc:`IndexError` is raised.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        allow_outside_chromosome : :py:class:`bool`
            If False, does not allow the interval to be outside the range of the chromosome.
            Attempting to pass an interval outside the chromosome range will raise an :py:exc:`IndexError`.
            If not specified, defaults to True. Any base outside the chromosome range will be returned as an 'N'.

        Returns
        -------
        :py:class:`str`
            The DNA sequence.
        """
        return mock_result(str)

    @mock
    @property
    def reference_genome(self):  # pragma: no cover
        """The name of the reference genome.

        Returns
        -------
        :py:class:`str`
           The name of the reference genome from which DNA will be extracted, e.g. "hg38".
        """
        return mock_result(str)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to the file from which DNA is extracted.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/data/hg19.2bit"
        """
        return mock_result(str)

    def __repr__(self):
        return "<GenomeDNA '{}'>".format(self.filename)
