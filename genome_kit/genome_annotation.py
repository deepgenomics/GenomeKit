# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import pickle
from typing import Optional

from . import _cxx
from ._cxx_util import mock, mock_result, mock_unreachable, strip_mock_bases
from .interval import Interval


#########################################################################


@strip_mock_bases
@_cxx.register
class Gene(_cxx.Gene, Interval):
    """An annotated gene.

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
        """The interval spanned by this gene.

        Note that :py:class:`~genome_kit.Gene` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this gene.
        """
        return mock_result(Interval)

    @mock
    @property
    def id(self):  # pragma: no cover
        """The ID of this gene.

        Returns either an Ensembl ID (for GENCODE) or a Entrez ID (for RefSeq)
        depending on which annotation is being represented.

        Returns
        -------
        :py:class:`str`
            The ID of this transcript, e.g. "ENSG00000115275.11" or "7841"
        """
        return mock_result(str)

    @mock
    @property
    def name(self):  # pragma: no cover
        """The name of this gene.

        Returns
        -------
        :py:class:`str`
            The name of this gene, e.g. "BRCA1"
        """
        return mock_result(str)

    @mock
    @property
    def type(self):  # pragma: no cover
        """The type of this gene.

        This represents the `gene_type` field in most GFF3 annotations.

        For GENCODE/Ensembl/NCBI RefSeq, the result can be any standard `GENCODE biotype
        <https://www.gencodegenes.org/pages/biotypes.html>`_

        For UCSC RefSeq, the result is "protein_coding" if any transcript is
        protein coding, and otherwise "non_coding".

        Returns
        -------
        :py:class:`str`
            The type of this gene, e.g. "protein_coding" or "non_coding"
        """
        return mock_result(str)

    @mock
    @property
    def level(self):  # pragma: no cover
        """The evidence level of this gene.

        This represents the `level` field in most GFF3 annotations,
        either 1, 2, or 3.

        Available for `gencode` and `ensembl` only.

        Returns
        -------
        :py:class:`int` | :py:data:`None`
            The evidence level of this gene.
        """
        return mock_result(int)

    @mock
    @property
    def transcripts(self):  # pragma: no cover
        """The transcripts for this gene.

        A shorthand property `trans` is also available.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            A copy of the list of annotated transcripts for this gene.
        """
        mock_unreachable()
        return [Transcript()]

    def __getstate__(self) -> bytes:
        genome = self.annotation_genome
        return pickle.dumps([genome, genome.genes.index_of(self)])

    def __setstate__(self, state: bytes) -> None:
        genome, index = pickle.loads(state)
        self._setstate(
            genome.genes, index
        )

    def __repr__(self):
        return "<Gene {} ({})>".format(self.id, self.name)


@_cxx.register
class GeneTable(_cxx.GeneTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all genes that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Gene()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all genes that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Gene()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all genes that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Gene()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all genes that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Gene()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all genes that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Gene()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all genes that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that overlap `interval`.
        """
        mock_unreachable()
        return [Gene()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all genes that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [gene for gene in genes if gene.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Gene`
            All genes that span exactly `interval`.
        """
        mock_unreachable()
        return [Gene()]

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
        """Lookup a gene by index or ID string.

        If `index` is a string the gene with matching ``id`` is returned (by linear search)::

            gene = genes['ENSG00000115275']  # MOGS gene

        If the ID string has a version suffix (``'ENSG00000115275.11'``) then the match must be exact.

        One can also iterate over all genes::

            for gene in genes:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int` | :py:class:`str`
            The integer index or ID string of the requested gene.

        Returns
        -------
        :py:class:`~genome_kit.Gene`
           The gene object identified by the given index.
        """
        return mock_result(Gene)

    def first_by_name(self, name: str) -> Optional[Gene]:
        """
        Returns
        -------
        :py:class:`~genome_kit.Gene`
            The first gene with the specified name or None if none matched.
        """
        return next((gene for gene in self if gene.name == name), None)

    def __repr__(self):
        return "<GeneTable, len() = {}>".format(len(self))


########################################################################


@strip_mock_bases
@_cxx.register
class Transcript(_cxx.Tran, Interval):
    """An annotated transcript.

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
        """The interval spanned by this transcript.

        Note that :py:class:`~genome_kit.Transcript` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this transcript.
        """
        return mock_result(Interval)

    @mock
    @property
    def id(self):  # pragma: no cover
        """The ID of this transcript.

        Returns either an Ensembl ID or a RefSeq ID depending on which
        annotation is being represented.

        Returns
        -------
        :py:class:`str`
            The ID of this transcript, e.g. "ENST00000233616.8" or "NM_006302"
        """
        return mock_result(str)

    @mock
    @property
    def tsl(self):  # pragma: no cover
        """The support level of this transcript.

        * tsl1: all splice junctions supported by at least one non-suspect mRNA
        * tsl2: best supporting mRNA is flagged as suspect or support from multiple ESTs
        * tsl3: only support is from a single EST
        * tsl4: best supporting EST is flagged as suspect
        * tsl5: no single transcript supports the model structure

        See `Ensemble transcript support level <http://www.ensembl.org/Help/Glossary?id=492>`_.

        Available for `gencode` only.

        Returns
        -------
        :py:class:`int` | :py:data:`None`
            The support level of this transcript, where values 1..5 and 6 correspond to
            GENCODE identifiers [tsl1, tsl2, tsl3, tsl4, tsl5, tslNA]
        """
        return mock_result(int)

    @mock
    @property
    def type(self):  # pragma: no cover
        """The type of this transcript.

        This represents the `transcript_type` field in most GFF3 annotations,
        which may differ from the parent gene's `type`.

        For GENCODE/Ensembl/NCBI RefSeq, the result can be any standard `GENCODE biotype
        <https://www.gencodegenes.org/pages/biotypes.html>`_

        For UCSC RefSeq, the result is either "protein_coding" or "non_coding".

        Returns
        -------
        :py:class:`str`
            The type of this transcript, e.g. "protein_coding" or "non_coding"
        """
        return mock_result(str)

    @mock
    @property
    def level(self):  # pragma: no cover
        """The evidence level of this transcript.

        This represents the `level` field in most GFF3 annotations,
        either 1, 2, or 3.

        Available for `gencode` and `ensembl` only.

        Returns
        -------
        :py:class:`int` | :py:data:`None`
            The evidence level of this transcript.
        """
        return mock_result(int)

    @mock
    @property
    def ccds_id(self):  # pragma: no cover
        """The CCDSID of this transcript's coding sequence, if applicable.

        Available for `gencode` only.

        Returns
        -------
        :py:class:`str` | :py:data:`None`
            The CCDSID of this transcript, e.g. "CCDS53759.1"
        """
        return mock_result(str)

    @mock
    @property
    def protein_id(self):  # pragma: no cover
        """The protein ID of this transcript, if applicable.

        Available for `gencode`, `ensembl`, `ucsc_refseq`, and `ncbi_refseq`.

        Returns
        -------
        :py:class:`str` | :py:data:`None`
            The protein ID of this transcript, e.g. "ENSP00000233616.4" or "NP_006293"
        """
        return mock_result(str)

    @mock
    @property
    def product(self):  # pragma: no cover
        """A description of the transcript's product.

        Available for `ucsc_refseq` and `ncbi_refseq` only.

        Returns
        -------
        :py:class:`str` | :py:data:`None`
            The product of this transcript, e.g. "breast cancer type 1 susceptibility protein isoform 2"
        """
        return mock_result(str)

    @mock
    @property
    def gene(self):  # pragma: no cover
        """The parent gene of this transcript.

        Returns
        -------
        :py:class:`~genome_kit.Gene`
            The parent gene for this transcript.
        """
        return mock_result(Gene)

    @mock
    @property
    def exons(self):  # pragma: no cover
        """The exons for this transcript.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            A copy of the list of annotated exons for this transcript.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    @property
    def introns(self):  # pragma: no cover
        """The introns for this transcript.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            A copy of the list of annotated introns for this transcript.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    @property
    def cdss(self):  # pragma: no cover
        """The coding sequences (CDSs) for this transcript.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            A copy of the list of annotated coding sequences for this transcript.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    @property
    def utr5s(self):  # pragma: no cover
        """The UTR5 segments for this transcript.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Utr`
            A copy of the list of annotated UTR5 segments for this transcript.
        """
        mock_unreachable()
        return [Utr()]

    @mock
    @property
    def utr3s(self):  # pragma: no cover
        """The UTR3 segments for this transcript.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Utr`
            A copy of the list of annotated UTR3 segments for this transcript.
        """
        mock_unreachable()
        return [Utr()]

    def __getstate__(self) -> bytes:
        genome = self.annotation_genome
        return pickle.dumps([genome, genome.transcripts.index_of(self)])

    def __setstate__(self, state: bytes) -> None:
        genome, index = pickle.loads(state)
        self._setstate(
            genome.transcripts, index
        )


    def __repr__(self):
        return "<Transcript {} of {}>".format(self.id, self.gene.name)


@_cxx.register
class TranscriptTable(_cxx.TranTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all transcripts that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Transcript()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all transcripts that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Transcript()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all transcripts that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Transcript()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all transcripts that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Transcript()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all transcripts that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Transcript()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all transcripts that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that overlap `interval`.
        """
        mock_unreachable()
        return [Transcript()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all transcripts that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [tran for tran in transcripts if tran.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            All transcripts that span exactly `interval`.
        """
        mock_unreachable()
        return [Transcript()]

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
        """Lookup a transcript by index or ID string.

        If `index` is a string the transcript with matching ``id`` is returned (by linear search)::

            tran = transcripts['ENST00000233616']  # MOGS transcript

        If the ID string has a version suffix (``'ENST00000233616.8'``) then the match must be exact.

        One can also iterate over all transcripts::

            for tran in transcripts:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int` | :py:class:`str`
            The integer index or ID string of the requested transcript.

        Returns
        -------
        :py:class:`~genome_kit.Transcript`
           The transcript object identified by the given index.
        """
        return mock_result(Transcript)

    def __repr__(self):
        return "<TranscriptTable, len = {}>".format(len(self))


########################################################################


@strip_mock_bases
@_cxx.register
class Exon(_cxx.Exon, Interval):
    """An annotated exon.

    Bases: :py:class:`~genome_kit.Interval`

    Each `Exon` object is equal only to itself. To build a set of unique exonic
    intervals, use the :py:attr:`~genome_kit.Exon.interval` attribute::

        >>> from genome_kit import Genome
        >>> genome = Genome('gencode.v19')
        >>> exonic_intervals = set([exon.interval for exon in genome.exons])

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
        """The interval spanned by this exon.

        Note that :py:class:`~genome_kit.Exon` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this exon.
        """
        return mock_result(Interval)

    @mock
    @property
    def id(self):  # pragma: no cover
        """The ID of this exon.

        Available for `gencode` and `ensembl` only.

        Returns
        -------
        :py:class:`str` | :py:data:`None`
            The ID of this exon, e.g. "ENSE00000846737.3"
        """
        return mock_result(str)

    @mock
    @property
    def index(self):  # pragma: no cover
        """The index of this exon within the transcript.

        Returns
        -------
        :py:class:`int`
            The index of the exon number within the parent transcript (0-based).
        """
        return mock_result(int)

    @mock
    @property
    def transcript(self):  # pragma: no cover
        """The parent transcript of this exon.

        A shorthand property `tran` is also available.

        Returns
        -------
        :py:class:`~genome_kit.Transcript`
            The parent transcript of this exon.
        """
        return mock_result(Gene)

    @mock
    @property
    def cds(self):  # pragma: no cover
        """The CDS on this exon, if any.

        Returns
        -------
        :py:class:`~genome_kit.Cds` | :py:data:`None`
            The CDS contained within this exon, if any.
        """
        return mock_result(Cds)

    @mock
    @property
    def utr5(self):  # pragma: no cover
        """The UTR5 on this exon, if any.

        Returns
        -------
        :py:class:`~genome_kit.Utr5` | :py:data:`None`
            The UTR5 contained within this exon, if any.
        """
        return mock_result(Utr)

    @mock
    @property
    def utr3(self):  # pragma: no cover
        """The UTR3 on this exon, if any.

        Returns
        -------
        :py:class:`~genome_kit.Utr3` | :py:data:`None`
            The UTR3 contained within this exon, if any.
        """
        return mock_result(Utr)

    @mock
    @property
    def next_exon(self):  # pragma: no cover
        """The next exon in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Exon` | :py:data:`None`
            The next annotated exon on the parent transcript, if any.
        """
        return mock_result(Exon)

    @mock
    @property
    def prev_exon(self):  # pragma: no cover
        """The previous exon in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Exon` | :py:data:`None`
            The previous annotated exon on the parent transcript, if any.
        """
        return mock_result(Exon)

    @mock
    @property
    def prev_intron(self):  # pragma: no cover
        """Returns the previous intron in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Intron` | :py:data:`None`
            The previous annotated intron on the parent transcript, if any.
        """
        return mock_result(Intron)

    @mock
    @property
    def next_intron(self):  # pragma: no cover
        """Returns the next intron in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Intron` | :py:data:`None`
            The next annotated intron on the parent transcript, if any.
        """
        return mock_result(Intron)

    def __getstate__(self) -> bytes:
        genome = self.annotation_genome
        return pickle.dumps([genome, genome.exons.index_of(self)])

    def __setstate__(self, state: bytes) -> None:
        genome, index = pickle.loads(state)
        self._setstate(
            genome.exons, index
        )

    def __repr__(self):
        return "<Exon {}/{} of {}>".format(
            self.index + 1, len(self.transcript.exons), self.transcript.id
        )


@_cxx.register
class ExonTable(_cxx.ExonTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all exons that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all exons that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all exons that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all exons that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all exons that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all exons that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that overlap `interval`.
        """
        mock_unreachable()
        return [Exon()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all exons that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [exon for exon in exons if exon.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Exon`
            All exons that span exactly `interval`.
        """
        mock_unreachable()
        return [Exon()]

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
        """Lookup an exon by index.

        Allows iteration over all exons::

            for exon in exons:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested exon.

        Returns
        -------
        :py:class:`~genome_kit.Exon`
           The exon object identified by the given index.
        """
        return mock_result(Exon)

    def __repr__(self):
        return "<ExonTable, len = {}>".format(len(self))


########################################################################


@strip_mock_bases
@_cxx.register
class Intron(_cxx.Intr, Interval):
    """An annotated intron.

    Bases: :py:class:`~genome_kit.Interval`

    Each `Intron` object is equal only to itself. To build a set of unique intronic
    intervals, use the :py:attr:`~genome_kit.Intron.interval` attribute::

        >>> from genome_kit import Genome
        >>> genome = Genome('gencode.v19')
        >>> intronic_intervals = set([intron.interval for intron in genome.introns])

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
        """The interval spanned by this intron.

        Note that :py:class:`~genome_kit.Intron` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this intron.
        """
        return mock_result(Interval)

    @mock
    @property
    def index(self):  # pragma: no cover
        """The index of the intron within the transcript.

        Returns
        -------
        :py:class:`int`
            The index of the intron within the parent transcript (0-based).
        """
        return mock_result(int)

    @mock
    @property
    def transcript(self):  # pragma: no cover
        """The parent transcript of this intron.

        A shorthand property `tran` is also available.

        Returns
        -------
        :py:class:`~genome_kit.Transcript`
            The parent transcript of this intron.
        """
        return mock_result(Transcript)

    @mock
    @property
    def next_intron(self):  # pragma: no cover
        """The next intron in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Intron` | :py:data:`None`
            The next annotated intron on the parent transcript, if any.
        """
        return mock_result(Intron)

    @mock
    @property
    def prev_intron(self):  # pragma: no cover
        """The previous intron in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Intron` | :py:data:`None`
            The previous annotated intron on the parent transcript, if any.
        """
        return mock_result(Intron)

    @mock
    @property
    def prev_exon(self):  # pragma: no cover
        """Returns the previous exon in the transcript.

        Returns
        -------
        :py:class:`~genome_kit.Exon`
            The previous annotated exon on the parent transcript.
        """
        return mock_result(Exon)

    @mock
    @property
    def next_exon(self):  # pragma: no cover
        """Returns the next exon in the transcript.

        Returns
        -------
        :py:class:`~genome_kit.Exon`
            The next annotated exon on the parent transcript.
        """
        return mock_result(Exon)

    def __getstate__(self) -> bytes:
        genome = self.annotation_genome
        return pickle.dumps([genome, genome.introns.index_of(self)])

    def __setstate__(self, state: bytes) -> None:
        genome, index = pickle.loads(state)
        self._setstate(
            genome.introns, index
        )

    def __repr__(self):
        return "<Intron {}/{} of {}>".format(
            self.index + 1, len(self.transcript.introns), self.transcript.id
        )


@_cxx.register
class IntronTable(_cxx.IntrTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all introns that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all introns that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all introns that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all introns that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all introns that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all introns that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that overlap `interval`.
        """
        mock_unreachable()
        return [Intron()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all introns that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [intr for intr in introns if intr.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Intron`
            All introns that span exactly `interval`.
        """
        mock_unreachable()
        return [Intron()]

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
        """Lookup an intron by index.

        Allows iteration over all introns::

            for intron in introns:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested intron.

        Returns
        -------
        :py:class:`~genome_kit.Intron`
           The intron object identified by the given index.
        """
        return mock_result(Intron)

    def __repr__(self):
        return "<IntronTable, len = {}>".format(len(self))


########################################################################


@strip_mock_bases
@_cxx.register
class Cds(_cxx.Cds, Interval):
    """An annotated coding sequence (CDS) element within an exon.

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
        """The interval spanned by this CDS.

        Note that :py:class:`~genome_kit.Cds` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this CDS.
        """
        return mock_result(Interval)

    @mock
    @property
    def phase(self):  # pragma: no cover
        """The phase of this CDS's start position.

        The `phase` is the number of bases that should be trimmed from the 5' end of the
        CDS element to reach the start of next codon:

        - phase 0 means the 5p end is in-frame
        - phase 1 means the 5p end is 1 bp away from being in-frame.
        - phase 2 means the 5p end is 2 bp away from being in-frame.

        See the `GFF3 documentation <http://gmod.org/wiki/GFF3>` for a precise explanation.

        The `frame` attribute in UCSC RefSeq is not available.

        Returns
        -------
        :py:class:`int`
            The phase this CDS, either 0, 1, or 2.
        """
        return mock_result(int)

    @mock
    @property
    def exon(self):  # pragma: no cover
        """The parent exon of this CDS.

        Returns
        -------
        :py:class:`~genome_kit.Exon`
            The parent exon of this CDS.
        """
        return mock_result(Exon)

    @mock
    @property
    def transcript(self):  # pragma: no cover
        """The parent transcript of this CDS.

        A shorthand property `tran` is also available.

        Returns
        -------
        :py:class:`~genome_kit.Transcript`
            The parent transcript of this CDS.
        """
        return mock_result(Transcript)

    @mock
    @property
    def next_cds(self):  # pragma: no cover
        """The next CDS in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Cds` | :py:data:`None`
            The next annotated CDS on the parent transcript, if any.
        """
        return mock_result(Cds)

    @mock
    @property
    def prev_cds(self):  # pragma: no cover
        """The previous CDS in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Cds` | :py:data:`None`
            The previous annotated CDS on the parent transcript, if any.
        """
        return mock_result(Cds)

    def __getstate__(self) -> bytes:
        genome = self.annotation_genome
        return pickle.dumps([genome, genome.cdss.index_of(self)])

    def __setstate__(self, state: bytes) -> None:
        genome, index = pickle.loads(state)
        self._setstate(
            genome.cdss, index
        )

    def __repr__(self):
        return "<Cds in Exon {}/{} of {}>".format(
            self.exon.index + 1, len(self.transcript.exons), self.transcript.id
        )


@_cxx.register
class CdsTable(_cxx.CdsTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    def find_5p_aligned(self, interval):  # pragma: no cover
        """Returns all CDSs that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that have 5' end aligned to the 5' end of `interval`.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    def find_3p_aligned(self, interval):  # pragma: no cover
        """Returns all CDSs that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that have 3' end aligned to the 3' end of `interval`.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    def find_5p_within(self, interval):  # pragma: no cover
        """Returns all CDSs that have 5'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds.end5.expand(0, 1) in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that have 5'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    def find_3p_within(self, interval):  # pragma: no cover
        """Returns all CDSs that have 3'-most position within `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds.end3.expand(1, 0) in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that have 3'-most position falling entirely within `interval`.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    def find_within(self, interval):  # pragma: no cover
        """Returns all CDSs that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that fall entirely within `interval`.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    def find_overlapping(self, interval):  # pragma: no cover
        """Returns all CDSs that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that overlap `interval`.
        """
        mock_unreachable()
        return [Cds()]

    @mock
    def find_exact(self, interval):  # pragma: no cover
        """Returns all CDSs that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [cds for cds in cdss if cds.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Cds`
            All coding sequences (CDSs) that span exactly `interval`.
        """
        mock_unreachable()
        return [Cds()]

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
        """Lookup a coding sequence element by index.

        Allows iteration over all coding sequences::

            for cds in cdss:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested CDS.

        The results will be ordered as in the source file.

        Returns
        -------
        :py:class:`~genome_kit.Cds`
           The CDS identified by the given index.
        """
        return mock_result(Cds)

    def __repr__(self):
        return "<CdsTable, len = {}>".format(len(self))


########################################################################

@strip_mock_bases
@_cxx.register
class Utr(_cxx.Utr, Interval):
    """An annotated UTR element within an exon.

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
        """The interval spanned by this UTR segment.

        Note that :py:class:`~genome_kit.Utr` also inherits from :py:class:`~genome_kit.Interval`.

        Returns
        -------
        :py:class:`~genome_kit.Interval`
            The interval spanned by this UTR segment.
        """
        return mock_result(Interval)

    @mock
    @property
    def exon(self):  # pragma: no cover
        """The parent exon of this UTR segment.

        Returns
        -------
        :py:class:`~genome_kit.Exon`
            The parent exon of this UTR segment.
        """
        return mock_result(Exon)

    @mock
    @property
    def transcript(self):  # pragma: no cover
        """The parent transcript of this UTR segment.

        A shorthand property `tran` is also available.

        Returns
        -------
        :py:class:`~genome_kit.Transcript`
            The parent transcript of this UTR segment.
        """
        return mock_result(Transcript)

    @mock
    @property
    def next_utr(self):  # pragma: no cover
        """The next UTR in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Utr` | :py:data:`None`
            The next annotated UTR segment on the parent transcript, if any.
        """
        return mock_result(Utr)

    @mock
    @property
    def prev_utr(self):  # pragma: no cover
        """The previous UTR in the transcript, if any.

        Returns
        -------
        :py:class:`~genome_kit.Utr` | :py:data:`None`
            The previous annotated UTR segment on the parent transcript, if any.
        """
        return mock_result(Utr)

    _FIVE_PRIME = "5"
    _THREE_PRIME = "3"

    def __getstate__(self) -> bytes:
        genome = self.annotation_genome
        try:
            return pickle.dumps([genome, genome.utr5s.index_of(self), Utr._FIVE_PRIME])
        except ValueError:
            return pickle.dumps([genome, genome.utr3s.index_of(self), Utr._THREE_PRIME])

    def __setstate__(self, state: bytes) -> None:
        genome, index, side = pickle.loads(state)
        if side == Utr._FIVE_PRIME:
            self._setstate(
                genome.utr5s, index
            )
        else:
            self._setstate(
                genome.utr3s, index
            )

    def __repr__(self):
        return "<Utr in Exon {}/{} of {}>".format(
            self.exon.index + 1, len(self.transcript.exons), self.transcript.id
        )


@_cxx.register
class UtrTable(_cxx.UtrTable):
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

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
        """Lookup a UTR element by index.

        Allows iteration over all coding sequences::

            for utr in utr5s:
                ...

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested UTR.

        The results will be ordered as in the source file.

        Returns
        -------
        :py:class:`~genome_kit.Utr`
           The UTR segment identified by the given index.
        """
        return mock_result(Utr)

    def __repr__(self):
        return "<UtrTable, len = {}>".format(len(self))

########################################################################

@_cxx.register
class GenomeAnnotation(_cxx.GenomeAnno):
    """Annotations for a reference genome.

    This object should not be directly created.
    Instead, create a :py:class:`~Genome` object::

        genome = Genome("gencode.v19")
        for gene in genome.annotation.genes:
            ...

    The :py:class:`~genome_kit.Genome` object also exposes shortcut
    attributes, for convenience::

        for gene in genome.genes:
            ...

    See :py:class:`~genome_kit.Genome` for a list of annotations
    available, e.g. ``gencode.v29``, ``ucsc_refseq.2017-06-25``, etc.

    .. note:: **UCSC RefSeq can contain multiple versions of the same gene/transcript.**
       These tend to be in ambiguously mapped regions. In GenomeKit, they're
       retained as separate entries in the gene and transcript tables, with
       distinct intervals but otherwise identical.

       For example, in GENCODE v19 the SMN1/2 genes are located on chromosome
       5 at 70220767-70249769 and 69345349-69374349 respectively.
       In UCSC RefSeq, SMN1 is duplicated at both 70220767-70248838 and
       69345349-69373418, so there are two gene entries named SMN1, each
       with their own version of transcripts NM_000344, NM_001297715, NM_022874.
       Similarly, SMN2 is duplicated at both 70220767-70248842 and 69345349-69373422.
    """

    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    @mock
    @property
    def genes(self):  # pragma: no cover
        """An indexed table of annotated genes.

        Returns
        -------
        :py:class:`~genome_kit.GeneTable`
           An object that supports looping or queries over annotated genes.
        """
        return mock_result(GeneTable)

    @mock
    @property
    def transcripts(self):  # pragma: no cover
        """An indexed table of annotated transcripts.

        A shorthand property `trans` is also available.

        Returns
        -------
        :py:class:`~genome_kit.TranscriptTable`
           An object that supports looping or queries over annotated transcripts.
        """
        return mock_result(TranscriptTable)

    @mock
    @property
    def exons(self):  # pragma: no cover
        """An indexed table of annotated exons.

        Returns
        -------
        :py:class:`~genome_kit.ExonTable`
           An object that supports looping or queries over annotated exons.
        """
        return mock_result(ExonTable)

    @mock
    @property
    def introns(self):  # pragma: no cover
        """An indexed table of annotated introns.

        Returns
        -------
        :py:class:`~genome_kit.IntronTable`
           An object that supports looping or queries over annotated introns.
        """
        return mock_result(IntronTable)

    @mock
    @property
    def cdss(self):  # pragma: no cover
        """An indexed table of annotated coding sequences (CDSs) within exons.

        This is shorthand for `annotation.cdss`.

        Returns
        -------
        :py:class:`~genome_kit.CdsTable`
           An object that supports looping or queries over annotated coding sequences.
        """
        return mock_result(CdsTable)

    @mock
    @property
    def utr5s(self):  # pragma: no cover
        """An indexed table of UTR5s within exons.

        This is shorthand for `annotation.utr5s`.

        Returns
        -------
        :py:class:`~genome_kit.UtrTable`
           An object that supports looping or queries over annotated coding sequences.
        """
        return mock_result(UtrTable)

    @mock
    @property
    def utr3s(self):  # pragma: no cover
        """An indexed table of UTR3s within exons.

        This is shorthand for `annotation.utr3s`.

        Returns
        -------
        :py:class:`~genome_kit.UtrTable`
           An object that supports looping or queries over annotated coding sequences.
        """
        return mock_result(UtrTable)

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to the file from which annotations are retrieved.

        Returns
        -------
        :py:class:`str`
           The path to the file, e.g. "/data/gencode.v19.annotation.dganno"
        """
        return mock_result(str)

    @mock
    @staticmethod
    def binary_version():  # pragma: no cover
        """The binary file format version that this build supports.

        Returns
        -------
        :py:class:`int`
           The version number.
        """
        return mock_result(int)

    @mock
    @staticmethod
    def build_gencode(infile, outfile, reference_genome):  # pragma: no cover
        """Build a binary version of GENCODE/Ensembl annotations from a GFF3 file.

        Currently supports conversion from `gff3[.gz]`.

        In addition to the versioned `dganno` file, a `cfg` file will also be created.

        Parameters
        ----------
        infile : :py:class:`str`
            The path to the original GFF3 file.

        outfile : :py:class:`str`
            The path to the destination `dganno` file.

        reference_genome : :class:`~genome_kit.Genome`
            The reference assembly used to annotate ``infile``.

        Returns
        -------
        :class:`list` of :class:`str`
            A list of built files created from the GFF3.
        """
        return mock_result(list)

    @mock
    @staticmethod
    def build_ncbi_refseq(infile, outfile, reference_genome):  # pragma: no cover
        """Build a binary version of NCBI annotations from a GFF3 file.

        Currently supports conversion from `gff3[.gz]`.

        In addition to the versioned `dganno` file, a `cfg` file will also be created.

        Parameters
        ----------
        infile : :py:class:`str`
            The path to the original GFF3 file.

        outfile : :py:class:`str`
            The path to the destination `dganno` file.

        reference_genome : :class:`~genome_kit.Genome`
            The reference assembly used to annotate ``infile``.

        Returns
        -------
        :class:`list` of :class:`str`
            A list of built files created from the GFF3.
        """
        return mock_result(list)

    @mock
    @staticmethod
    def build_ucsc_refseq(ucsc_db_dir, outfile, reference_genome):  # pragma: no cover
        """Build a binary version of UCSC RefSeq annotations.

        The `ucsc_db_dir` should be a directory that already contains
        the following files:

        * ``refGene.txt[.gz]``
        * ``refLink.txt[.gz]``

        These files will be parsed and used to build a `dganno` file.

        In addition to the versioned `dganno` file, a `cfg` file will also be created.

        Parameters
        ----------
        ucsc_db_dir : :py:class:`str`
            The path to a directory with the UCSC RefSeq source files.

        outfile : :py:class:`str`
            The path to the destination `dganno` file.

        reference_genome : :class:`~genome_kit.Genome`
            The reference assembly used  ``ucsc_db_dir``.

        Returns
        -------
        :class:`list` of :class:`str`
            A list of built files created from the UCSC RefSeq files.
        """
        return mock_result(list)

    def __repr__(self):
        return '<GenomeAnnotation "{}">'.format(self.filename)


########################################################################
