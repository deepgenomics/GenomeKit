# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import os
from weakref import WeakValueDictionary

from . import _cxx
from . import genome_annotation as _ga
from . import interval as _interval
from ._apply_variants import apply_variants, check_variants_list
from ._cxx_util import mock, mock_result
from ._util import reverse_complement
from .variant import Variant


@_cxx.register
class Genome(_cxx.Genome):
    """Initialize a reference genome and its corresponding resources (DNA, annotations, etc).

    A genome object represents a reference genome and resources available for it,
    such as DNA and annotations::

        >>> genome = Genome("gencode.v19")
        >>> exon = genome.exons[0]                                # 1st annotated exon.
        >>> dna  = genome.dna(exon)                               # Extract DNA.

    The initializer argument specifies a versioned resource::

        # Access resources that are hg19
        Genome("hg19")

        # GENCODE v29 implies hg38
        Genome("gencode.v29")

    The currently recognized strings depend on the data made available by data_manager.
    Some common ones are::

        reference_genomes = [
            "hg19",
            "hg19.p13.plusMT",
            "hg38",
            "hg38.p12",
            "hg38.p13",
            "hg38.p14",
            "mm10.p6",
            "mm39",
            "rn6",
            "macFas5",
            "susScr11",
        ]

        annotations = [
            "gencode.v19",                  # implies hg19.p13.plusMT
            "gencode.v29",                  # implies hg38.p12
            "gencode.v29.basic",            # implies hg38.p12
            "gencode.v29lift37",            # implies hg19
            "gencode.v29lift37.basic",      # implies hg19
            "gencode.v41",                  # implies hg38.p13
            "gencode.v41.basic",            # implies hg38.p13
            "gencode.v41lift37",            # implies hg19
            "gencode.v41lift37.basic",      # implies hg19
            "ucsc_refseq.2017-06-25",       # implies hg19
            "ncbi_refseq.v105.20190906",    # implies hg19.p13.plusMT
            "ncbi_refseq.v109",             # implies hg38.p12
            "ncbi_refseq.v110",             # implies hg38.p14
            "gencode.vM19",                 # implies mm10.p6
            "gencode.vM19.basic",           # implies mm10.p6
            "gencode.vM30",                 # implies mm39
            "gencode.vM30.basic",           # implies mm39
            "gencode.vM31",                 # implies mm39
            "gencode.vM31.basic",           # implies mm39
            "ncbi_refseq.m39.v109",         # implies mm39
            "ensembl.Rnor_6.0.88",          # implies rn6
            "ncbi_refseq.Macfas_5.0.v101",  # implies macFas5
            "ensembl.Macfas_5.0.95",        # implies macFas5
            "ensembl.Sscrofa11.1.98",       # implies susScr11
        ]

    If a resource is requested after object creation, but the version
    needed cannot be disambiguated from `config`, an exception is raised.

    Parameters
    ----------
    config : :py:class:`str`
        Identifies a versioned resource that this genome object is
        expected to serve.
    """

    def __init__(self, config):  # pragma: no cover
        self._chromosome_sizes = None
        self._appris_indices = None
        self._appris_transcripts = None
        self._appris_transcripts_by_gene = None
        self._appris_principality_strings = ('PRINCIPAL:1', 'PRINCIPAL:2', 'PRINCIPAL:3', 'PRINCIPAL:4', 'PRINCIPAL:5',
                                             'ALTERNATIVE:1', 'ALTERNATIVE:2')

    genome_cache = WeakValueDictionary()

    def __reduce__(self):
        return (type(self), (self.config, ))

    def __new__(cls, config):  # pragma: no cover
        if cls == Genome:
            genome = Genome.genome_cache.get(config, None)
            if genome:
                return genome

        genome = _cxx.Genome.__new__(cls, config)

        if cls == Genome:  # do not cache subclasses (like mocks)
            Genome.genome_cache[config] = genome

        return genome

    def interval(self, chromosome, strand, start, end):
        """Shorthand method for creating an Interval with the same refg as the genome object called on.
        """
        return _interval.Interval(chromosome, strand, start, end, self.refg)

    def variant(self, variant_string):
        """Shorthand method for creating a Variant with the same refg as the genome object called on.
        """
        return Variant.from_string(variant_string, self)

    @mock
    def dna(self, interval, allow_outside_chromosome=True):  # pragma: no cover
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

    def variant_dna(self, interval, variants):
        """Extract DNA from a genome with specific variants applied to it.

        Follows the same rules as :py:meth:`~genome_kit.VariantGenome.dna`.
        Intended for applying sets of variants that may overlap `interval`.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        variants : :py:class:`~genome_kit.Variant` | :py:class:`list` of (:py:class:`~genome_kit.Variant`)
            Either a single variant or a list of variants. See :py:class:`~genome_kit.VariantGenome`.

        Returns
        -------
        :py:class:`str`
            The variant DNA sequence.
        """
        return apply_variants(self.dna, check_variants_list(self, variants), interval)

    @mock
    @property
    def annotation(self):  # pragma: no cover
        """Access to annotations.

        A shorthand property `anno` is also available.

        Returns
        -------
        :py:class:`~genome_kit.GenomeAnnotation`
            An object that provides access to annotated genes, transcripts, and exons.
        """
        return mock_result(_ga.GenomeAnnotation)

    @mock
    @property
    def genes(self):  # pragma: no cover
        """Access to annotated genes.

        This is shorthand for `annotation.genes`.

        Returns
        -------
        :py:class:`~genome_kit.GeneTable`
           An object that supports looping or queries over annotated genes.
        """
        return mock_result(_ga.GeneTable)

    @mock
    @property
    def transcripts(self):  # pragma: no cover
        """Access to annotated transcripts.

        This is shorthand for `annotation.transcripts`.

        A shorthand property `trans` is also available.

        Returns
        -------
        :py:class:`~genome_kit.TranscriptTable`
           An object that supports looping or queries over annotated transcripts.
        """
        return mock_result(_ga.TranscriptTable)

    @mock
    @property
    def exons(self):  # pragma: no cover
        """Access to annotated exons.

        This is shorthand for `annotation.exons`.

        Returns
        -------
        :py:class:`~genome_kit.ExonTable`
           An object that supports looping or queries over annotated exons.
        """
        return mock_result(_ga.ExonTable)

    @mock
    @property
    def introns(self):  # pragma: no cover
        """Access to annotated introns.

        This is shorthand for `annotation.introns`.

        Returns
        -------
        :py:class:`~genome_kit.IntronTable`
           An object that supports looping or queries over annotated introns.
        """
        return mock_result(_ga.IntronTable)

    @mock
    @property
    def cdss(self):  # pragma: no cover
        """Access to annotated coding sequences (CDSs).

        This is shorthand for `annotation.cdss`.

        Returns
        -------
        :py:class:`~genome_kit.CdsTable`
           An object that supports looping or queries over annotated coding sequences.
        """
        return mock_result(_ga.CdsTable)

    @mock
    @property
    def utr5s(self):  # pragma: no cover
        """Access to annotated UTR5 sequences.

        This is shorthand for `annotation.utr5s`.

        Returns
        -------
        :py:class:`~genome_kit.Utr5Table`
           An object that supports looping or queries over annotated coding sequences.
        """
        return mock_result(_ga.Utr5Table)

    @mock
    @property
    def utr3s(self):  # pragma: no cover
        """Access to annotated UTR3 sequences.

        This is shorthand for `annotation.utr3s`.

        Returns
        -------
        :py:class:`~genome_kit.Utr3Table`
           An object that supports looping or queries over annotated coding sequences.
        """
        return mock_result(_ga.Utr3Table)

    def find_motif(self, interval, motif, match_position=0, find_overlapping_motifs=False):
        """Find a genomic motif in an interval on a :py:class:`~genome_kit.Genome`.

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
        :py:class:`list` of :py:class:`~genome_kit.Interval`
            A list of all matches (length-0 intervals) found in the query interval.

        Example
        -------::
            >>> genome37 = Genome('hg19')
            >>> interval = genome37.interval('chr1', '+', 40000, 40080)
            >>> motif = 'AA'
            >>> genome37.find_motif(interval, motif)  # doctest:+NORMALIZE_WHITESPACE
            [Interval("chr1", "+", 40031, 40031, "hg19"),
             Interval("chr1", "+", 40061, 40061, "hg19")]
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
            motif = reverse_complement(motif)

        dna = self.dna(interval.as_positive_strand())

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

            pos_genome = dna_start + pos + match_position

            hit = _interval.Interval(interval.chromosome, interval.strand, pos_genome, pos_genome, interval.refg)
            motif_hits.append(hit)
            offset = pos + offset_step

        return motif_hits

    @mock
    @property
    def reference_genome(self):  # pragma: no cover
        """The name of the reference genome.

        Returns
        -------
        :py:class:`str`
           The name of the reference genome that any resources map to, e.g. "hg19" or "hg38".
        """
        return mock_result(str)

    @property
    def chromosomes(self):
        """A list of valid chromosome names for this genome."""
        return list(self._chrom_sizes().keys())

    def chromosome_size(self, chromosome):
        """Returns the size of the given chromosome, in bp.

        Returns
        -------
        :py:class:`int`
            The size of the given chromosome, as number of basepairs.
        """
        try:
            return self._chrom_sizes()[chromosome]
        except KeyError:
            error_message = "No such chromosome: %s. \nSupported values: %s" % (chromosome, self._chrom_sizes().keys())
            raise ValueError(error_message)

    @mock
    @property
    def config(self):
        """The configuration strings of this genome.

        Returns
        -------
        :py:class:`str`
            The config used to initialize this object.
        """
        return mock_result(str)

    @mock
    @property
    def data_dir(self):  # pragma: no cover
        """The data_dir from which resources are opened.

        Returns
        -------
        :py:class:`str`
            The path to where resource files are being opened.
        """
        return mock_result(str)

    def __repr__(self):
        return f'Genome("{self.config}")'

    def _chrom_sizes(self):
        if self._chromosome_sizes:
            return self._chromosome_sizes

        self._chromosome_sizes = {}

        with open(self._chrom_sizes_filename()) as file:
            for line in file:
                (key, val) = line.split()
                self._chromosome_sizes[str(key)] = int(val)

        return self._chromosome_sizes

    def _chrom_sizes_filename(self):
        from . import gk_data

        # If this Genome object has a custom data_dir resource, we must respect it.
        datafile_path = os.path.join(self.data_dir, f"{self.refg}.chrom.sizes")

        # Before returning, ensure that the path resolves to an actual file
        # on the local file system, pulling from store if necessary.
        return gk_data.resolve_datafile_path(datafile_path)

    def _load_appris(self):
        """Parse the APPRIS text file corresponding to this genome's annotations."""
        import pickle

        from . import gk_data
        from ._gk_data_config import get_appris_filename

        # Figure out which annotation version we should match with APPRIS
        datafile_name = get_appris_filename(self.config)
        datafile_path = os.path.join(self.data_dir, datafile_name)
        try:
            datafile_path = gk_data.resolve_datafile_path(datafile_path)
        except Exception as e:
            raise ApprisNotAvailableError(f"APPRIS not available for {self.annotation.filename}") from e

        with open(datafile_path, 'rb') as fp:
            data = pickle.load(fp)
        # List of int, where item [i] is the principality index (0..6) for self.transcripts[i]
        self._appris_indices = data[0]
        # List of Transcript indices, ordered by APPRIS index
        self._appris_transcripts = data[1]
        # List of lists, where item [i] is a list of Transcript indices belonging to self.genes[i],
        # ordered by APPRIS index.
        self._appris_transcripts_by_gene = data[2]

    def appris_transcripts(self, gene=None):
        """Returns a list of transcript objects ordered by our custom APPRIS index.
        If a gene is provided, returns a list of transcripts for that gene ordered by decreasing priority.

        Parameters
        ----------
        gene : :py:class:`~genome_kit.Gene` | :py:class:`str`, optional
            Restrict to a specific gene. If using `str`, use the geneID.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.Transcript`
            A list of transcripts that have APPRIS scores available. If a gene is provided as input, a list
            of transcripts for that gene in the order from most to least principle; so the first transcript should
            be used as the canonical.

        Raises
        ------
        ApprisNotAvailableError:
            APPRIS data is not available for this genome's annotations.
        """

        if self._appris_transcripts is None:
            self._load_appris()

        if gene is None:
            # Use all transcript indices ordered by principality.
            transcript_indices = self._appris_transcripts
        else:
            if isinstance(gene, str):
                # Get transcript indices ordered by APPRIS principality.
                gene = self.genes[gene]
            gene_idx = self.genes.index_of(gene)
            transcript_indices = self._appris_transcripts_by_gene[gene_idx]

        return [self.transcripts[i] for i in transcript_indices]

    def appris_principality(self, transcript, as_string=False):
        """Returns the APPRIS principality index (isoform specific number) for this transcript.
        An APPRIS principality index is an integer from 0..6 where 0 is considered most prinicpal and 6 considered least
        principal.
        In the original APPRIS database, the 7 principality levels are reported using strings, in the following order::

            appris_index_names = [
                "PRINCIPAL:1",     # index 0
                "PRINCIPAL:2",     # index 1
                "PRINCIPAL:3",     # index 2
                "PRINCIPAL:4",     # index 3
                "PRINCIPAL:5",     # index 4
                "ALTERNATIVE:1",   # index 5
                "ALTERNATIVE:2",   # index 6
            ]

        See the `APPRIS database documentation <http://appris.bioinfo.cnio.es/#/help/database/report>`_ for details.

        The APPRIS version defaults to `2017_05.v22` most annotations; see _build_appris.REMOTE_FILES for the rest.

        Parameters
        ----------
        transcript : :py:class:`~genome_kit.Transcript` | :py:class:`str`
            A transcript to query.
        as_string : :py:class:`bool`, optional
            Toggle True if the APPRIS principality string (e.g. `"PRINCIPAL:1"`) is desired.

        Returns
        -------
        :py:class:`int` | :py:class:`str`
            The APPRIS principality index [0..6] of this transcript, or the string representation if
            `as_string=True`

        Raises
        ------
        ApprisNotAvailableError:
            APPRIS data is not available for this genome's annotations.
        """
        if self._appris_indices is None:
            self._load_appris()

        if isinstance(transcript, str):
            transcript = self.transcripts[transcript]

        # Index of transcript object
        tidx = self.transcripts.index_of(transcript)

        if as_string:
            return self._appris_principality_strings[self._appris_indices[tidx]]

        return self._appris_indices[tidx]

class ApprisNotAvailableError(RuntimeError):
    pass