.. _quickstart:

----------
Quickstart
----------

Here is a quick tutorial covering GenomeKit's most important functionality.

Installation
------------

Install via conda
~~~~~~~~~~~~~~~~~

The best way to install GenomeKit is via the pre-compiled conda packages
available from https://anaconda.org/conda-forge/genomekit.

You can install GenomeKit with::

    mamba install genomekit


Basics
------

Initial Data Access
~~~~~~~~~~~~~~~~~~~

Data, such as assemblies, annotations, tracks, etc, are stored in custom-built binary files.

The files are required most commonly when creating :py:class:`~genome_kit.Genome` objects.
APIs for building these files are provided as part of the API. A selection of pre-built
data files is provided in a public Google Cloud Storage bucket, which is set as the default
data source.

The bucket is configured with `"requester pays" <https://cloud.google.com/storage/docs/requester-pays>`__
enabled, so you will need to `set up your GCloud credentials <https://cloud.google.com/sdk/docs/authorizing>`__
beforehand and set the ``GENOMEKIT_GCS_BILLING_PROJECT`` env var to your Google Cloud project.

For more details on creating your own data source, see `the section on Sharing data files <sharing data>`_.



Import the :py:mod:`genome_kit` package into your project::

    >>> import genome_kit as gk

For brevity, this tutorial assumes you've imported several core types
directly into your namespace::

    >>> from genome_kit import Genome
    >>> from genome_kit import Interval
    ...

Most GenomeKit code will need two core object types:

- :py:class:`~genome_kit.Interval` - an interval on a reference genome
- :py:class:`~genome_kit.Genome` - resources available for a reference genome

An `interval` is initialized from five values:
`chromosome`, `strand`, `start`, `end`, and `reference_genome`

    >>> interval = Interval("chr7", "+", 117120016, 117120201, "hg19")

Irrespective of strand, an interval always has `start <= end`.
The span of an interval excludes the `end` position, so the number of bases
spanned by an interval is `end -- start`.

    >>> len(interval)
    185

GenomeKit intervals can be described as "0-based, exclusive end",
whereas UCSC browser positions are "1-based, inclusive end"::

    >>> interval.as_ucsc()
    'chr7:117120017-117120201'

A `genome` provides convenient access to resources associated with a
reference genome::

    >>> genome = Genome("hg19")  # Equivalently "hg19"
    >>> genome.dna(interval)
    'AATTGGAAGCAAA...AACTTTTTTTCAG'

Some resources are versioned and are only accessible if the intended
version is unambiguous from the configuration strings. For example,
specifying a GENCODE version enables annotations::

    >>> genome = Genome("gencode.v19")             # Implies "hg19"
    >>> gene = genome.genes["ENSG00000001626.10"]  # Gene object
    >>> tran = gene.transcripts[2]                 # Transcript object
    >>> exon = tran.exons[0]                       # Exon object
    >>> genome.dna(exon)
    'AATTGGAAGCAAA...AACTTTTTTTCAG'

The above exon has the same interval that we defined earlier,
plus several other attributes::

    >>> exon.interval
    Interval("chr7", "+", 117120016, 117120201, "hg19")
    >>> exon.index
    0
    >>> exon.transcript
    <Transcript ENST00000003084.6 of CFTR>
    >>> exon.cds
    <Cds in Exon 1/27 of ENST00000003084.6>
    >>> exon.next_exon
    <Exon 2/27 of ENST00000003084.6>

Tracks
------

GenomeKit users can also build tracks via
:py:class:`~genome_kit.GenomeTrackBuilder`.

.. note::
    When building a track, data is ordered according to the strandedness argument
    passed to the builder's constructor:

    * ``"single_stranded"``: both strands share the same data. The data is applied in Interval coordinate (reference strand) order.
    * ``"strand_unaware"``: ignores the Interval strand, data is applied in Interval coordinate (reference strand) order.
    * ``"strand_aware"``: data is applied from 5" end to 3" end (sense strand order).

    >>> track = GenomeTrackBuilder("neg.gtrack", "u3", "strand_unaware", Genome("hg19"))
    >>> interval = Interval("chr1", "-", 10, 15, "hg38")
    >>> track.set_data(interval, np.arange(0, len(interval), dtype=np.uint8))
    >>> track.finalize()
    >>> track = GenomeTrack("neg.gtrack")
    >>> track(interval)
    array([[4],
           [3],
           [2],
           [1],
           [0]], dtype=uint8)
    >>> track = GenomeTrackBuilder("neg.gtrack", "u3", "strand_aware", Genome("hg19"))
    >>> track.set_data(interval, np.arange(0, len(interval), dtype=np.uint8))
    >>> track.finalize()
    >>> track = GenomeTrack("neg.gtrack")
    >>> track(interval)
    array([[0],
           [1],
           [2],
           [3],
           [4]], dtype=uint8)

Annotations
-----------

GenomeKit provides access to GENCODE and RefSeq annotations.

Walking the structure of an annotation is straight-forward::

    genome = Genome("gencode.v19")
    for gene in genome.genes:          # Each gene
        print(gene)
        for tran in gene.transcripts:  # Each transcript on the gene
            print("  ", tran)
            for exon in tran.exons:    # Each exon on the transcript
                print("     ", exon)

You can run the above example from the GenomeKit directory::

    $ python demos/walk_annotations.py
    <Gene ENSG00000223972.4 (DDX11L1)>
       <Transcript ENST00000456328.2 of DDX11L1>
          <Exon 1/3 of ENST00000456328.2>
          <Exon 2/3 of ENST00000456328.2>
          <Exon 3/3 of ENST00000456328.2>
    ...

Annotations are accessible via following attributes:
:py:attr:`~.genome_kit.Genome.genes`,
:py:attr:`~.genome_kit.Genome.transcripts`,
:py:attr:`~.genome_kit.Genome.exons`,
:py:attr:`~.genome_kit.Genome.introns`, and
:py:attr:`~.genome_kit.Genome.cdss`.
Each of these can be thought of as an indexed table, like in a database.
For example, `exons` is an instance of type
:py:class:`~.genome_kit.ExonTable` and each row in that table is an instance of
:py:class:`~.genome_kit.Exon`, with the fields you'd expect.

Most importantly, each annotation table is indexed for fast
positional queries::

    >>> # First exon of a CFTR transcript
    >>> exon = genome.transcripts["ENST00000003084.6"].exons[0]
    >>> genome.exons.find_overlapping(exon)
    [<Exon 1/27 of ENST00000003084.6>,
     <Exon 1/26 of ENST00000454343.1>,
     <Exon 1/26 of ENST00000426809.1>]

The following methods take a single `interval` argument and are available
for every element table:

- :py:meth:`~genome_kit.ExonTable.find_overlapping` - elements overlapping
  `interval`.

- :py:meth:`~genome_kit.ExonTable.find_within` - elements falling
  within `interval`.

- :py:meth:`~genome_kit.ExonTable.find_exact` - elements exactly spanning `interval`.

- :py:meth:`~genome_kit.ExonTable.find_5p_aligned` - elements with 5' end
  aligned to the 5' end of `interval`.

- :py:meth:`~genome_kit.ExonTable.find_3p_aligned` - elements with 3' end
  aligned to the 3' end of `interval`.

- :py:meth:`~genome_kit.ExonTable.find_5p_within` - elements with 5'-most
  position within `interval`.

- :py:meth:`~genome_kit.ExonTable.find_3p_within` - elements with 3'-most
  position within `interval`.

These methods are useful for mapping positions to annotation elements.
See :py:class:`~.genome_kit.GenomeAnnotation` for the top-level object that
ultimately owns all the annotation tables.




Working with Intervals
----------------------

GenomeKit uses the following convention for
:py:class:`~genome_kit.Interval`:

1. intervals are always stranded (+ or --),
2. positions are internally 0-based, and
3. the span of an interval excludes its `end` position.

Intervals are initialized using what is called the `DNA0` convention where,
irrespective of strand, where `start <= end`.
Intervals can be empty (`start == end`).
The resulting interval spans the same positions in the genome
as a Python slice operation `[start:end]`.

For example, an interval with `start=3` and `end=7` spans
the four over-lined positions below, irrespective of strand::

            [__________]
    0  1  2  3  4  5  6  7  8  9  10


To see what we can do with intervals, let's create a few::

     >>> #  0123456789
     >>> #  aaaaabbbbb
     >>> #     cccc
     >>> #     d
     >>> a = Interval("chr1", "+", 0,  5, "hg38")
     >>> b = Interval("chr1", "+", 5, 10, "hg38")
     >>> c = Interval("chr1", "+", 3,  7, "hg38")
     >>> d = Interval("chr1", "+", 3,  4, "hg38")

     >>> len(a), len(b), len(c), len(d)
     (5, 5, 4, 1)

     >>> a.contains(c), c.within(a), a.contains(d), d.within(a)
     (False, False, True, True)

     >>> a.overlaps(b), a.overlaps(c)
     (False, True)

     >>> a.upstream_of(b), b.dnstream_of(a)
     (True, True)
     >>> c.upstream_of(b), b.dnstream_of(c)
     (False, False)

     >>> a == b, a == d
     (False, False)
     >>> a != b, a != d
     (True, True)

Intervals on opposite strands effectively live in different universes::

    >>> #  0123456789
    >>> #  xxxxxyyyyy
    >>> #     zzzz
    >>> #     w
    >>> x = a.as_opposite_strand()
    >>> y = b.as_opposite_strand()
    >>> z = c.as_opposite_strand()
    >>> w = d.as_opposite_strand()

    >>> x
    Interval("chr1", "-", 0, 5, "hg38")
    >>> y
    Interval("chr1", "-", 5, 10, "hg38")
    >>> z
    Interval("chr1", "-", 3, 7, "hg38")
    >>> w
    Interval("chr1", "-", 3, 4, "hg38")

    >>> x.overlaps(d)  # opposite strands
    False
    >>> x.contains(d)  # opposite strands
    False
    >>> x.contains(w)
    True

    >>> x.upstream_of(y), y.upstream_of(x)
    (False, True)

    >>> x == a   # opposite strands
    False

Given an interval, you can also build new intervals centered around its
5' end (upstream) and 3' end (downstream),
which depends on strand::

    >>> interval = Interval("chr1", "-", 4,  8, "hg38")
    >>> interval.end5
    Interval("chr1", "-", 8, 8, "hg38")

    >>> interval.end3
    Interval("chr1", "-", 4, 4, "hg38")

    >>> interval.end3.expand(2, 3)
    Interval("chr1", "-", 1, 6, "hg38")

Notice that the `end5` and `end3` attributes return empty (length-0) intervals.
This is a general convention in GenomeKit: *an empty interval denotes the
space between two consecutive positions in the genome*.
This convention is useful for defining alignments or defining new intervals
relative to the empty one via :py:meth:`~genome_kit.Interval.expand`.
For example, the intervals in the above code example can be visualized
as follows::

               [__________]           # Interval("chr1", "-", 4,  8, "hg38")
    0  1  2  3  4  5  6  7  8  9  10

                          []          # interval.end5
    0  1  2  3  4  5  6  7  8  9  10

              []                      # interval.end3
    0  1  2  3  4  5  6  7  8  9  10

      [_____________]                 # interval.end3.expand(2, 3)
    0  1  2  3  4  5  6  7  8  9  10


Feature extraction basics
-------------------------

GenomeKit currently provides DNA sequence, accessible via the
`dna` attribute of :py:class:`.Genome`. The extracted DNA is
automatically reverse-complemented according to strand::

    >>> a = Interval("chr7", "+", 117120016, 117120201, "hg19")
    >>> b = a.as_opposite_strand()

    >>> genome = Genome("hg19")
    >>> genome.dna(a)
    'AATTGGAAGCAAA...AACTTTTTTTCAG'
    >>> genome.dna(b)
    'CTGAAAAAAAGTT...TTTGCTTCCAATT'

GenomeKit makes it easy to extract features that correspond to annotated
elements, and to filter those elements by their attributes.
For example, we could extract DNA from all acceptor sites that meet certain
criteria::

    def has_coding_seq(e): return e.cds is not None   # Has CDS?
    def has_good_level(e): return e.tran.level <= 2   # Is level 1 or 2?
    def not_first_exon(e): return e.index > 0         # Is not first exon?

    genome = Genome("gencode.v19")                    # 1196293 exons
    exons = filter(has_coding_seq, genome.exons)      #  724078 remaining
    exons = filter(has_good_level, exons)             #  605573 remaining
    exons = filter(not_first_exon, exons)             #  558099 remaining
    sites = set(exon.end5 for exon in exons)          #  198992 unique

    for site in sites:
        print(genome.dna(site.expand(5, 5)))

The above outputs the following 10nt sequences surrounding acceptor sites::

    TGCAGGGAAC   # Note they are all sense-strand (AG)
    TTCAGCTGCT   # because exon.end5 knows the strand.
    TGTAGGAAAC
    TCCAGGCTAT
    GCCAGAGGAC
    GACAGAACCA
    CCCAGATTGG
    ...

GenomeKit also make it easy to map positions to the nearest annotated
element. For example, we could map branch sites to their nearest downstream
acceptor site::

    genome = Genome("gencode.v19")

    # We want to map each of these to an acceptor at most 100nt downstream
    coords = [ Interval.from_dna0_coord("chr1", "-",  91661, "hg19"),
               Interval.from_dna0_coord("chr1", "-", 169295, "hg19"),
               Interval.from_dna0_coord("chr1", "+", 320861, "hg19") ]

    for coord in coords:
        window = coord.expand(0, 100)                # Find all exons with
        exons = genome.exons.find_5p_within(window)  # acceptor in window
        for exon in exons:
            print(coord.as_ucsc(), "-->", exon)

The above outputs the following mappings::

    chr1:91662-91662 --> <Exon 4/4 of ENST00000466430.1>
    chr1:169296-169296 --> <Exon 3/8 of ENST00000466557.2>
    chr1:320862-320862 --> <Exon 2/3 of ENST00000432964.1>  # 1st candidate
    chr1:320862-320862 --> <Exon 2/4 of ENST00000601486.1>  # 2nd candidate
    chr1:320862-320862 --> <Exon 1/3 of ENST00000599771.2>  # 3rd candidate

Variants
--------

Individual genomic variants are represented by a :py:class:`~.Variant`::

    from genome_kit import Variant

A variant is defined by a chromosome, 0-based position (DNA0), reference allele,
alternate allele, and reference genome::

    >>> variant = Variant("chr7", 117120148, "AT", "G", "hg19")
    >>> variant
    <Variant chr7:117120148:AT:G:hg19>

A variant is a subclass of :py:class:`~.Interval`, and also has an
:py:attr:`~.Variant.interval` attribute.
The interval spans the reference allele::

    >>> variant.start, variant.end, len(variant)
    (117120148, 117120150, 2)

    >>> variant.interval
    Interval("chr7", "+", 117120148, 117120150, "hg19")

Variants may also be created from a string where the position is 1-based (DNA1),
which is the convention of UCSC and Clinvar::

    >>> genome = Genome("hg19")
    >>> variant = genome.variant("chr7:117,120,149:AT:G")               # First way
    >>> variant = Variant.from_string('chr7:117,120,149:AT:G', genome)  # Second way
    >>> variant
    <Variant chr7:117120148:AT:G:hg19>

A variant created from a string is always validated against the reference genome,
raising an exception if the reference allele does not match.

ENSEMBL-style chromosome names (without leading ``chr``) are also allowed and are
automatically converted.

The GenomeKit variant format includes, but is more general than, the Clinvar variant convention when
it comes to insertion and deletions. Clinvar does not allow an empty ref or alt allele (requiring the
allele before the indel to be repeated), while empty alleles are allowed in GenomeKit. For example, the
variants ``chr7:117,120,150:A:-`` and ``chr7:117,120,149:CA:C`` are interpreted identically by GenomeKit.

Empty alleles can be specified by ``""``, ``"-"``, or ``"."``.


Variants from VCF Files
-----------------------

GenomeKit provides a binary VCF format that is compact, indexed by position.
The binary files can be opened by :py:class:`~genome_kit.VCFTable`, which returns
convenient :py:class:`~.Variant`-based objects::

    >>> from genome_kit import VCFTable

Suppose you have the following VCF saved as ``test.vcf.gz``::

    ##fileformat=VCFv4.2
    ##reference=GRCh37
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
    #CHROM POS     ID REF ALT QUAL FILTER INFO    FORMAT sample1  sample2  sample3
    1      949523  .  C   T   .    .      AF=0.00 GT:AD  0/0:0,1  0/1:0,2  0/0:0,3
    1      949608  .  G   A   .    .      AF=0.01 GT:AD  0/0:0,4  0/1:0,5  0/0:0,6
    1      949696  .  -   G   .    .      AF=0.02 GT:AD  0/0:0,7  0/1:0,8  0/1:0,9
    1      949739  .  G   TC  .    .      AF=0.03 GT:AD  0/1:0,10 0/0:0,11 1/1:0,12
    1      977028  .  G   T   .    .      AF=0.04 GT:AD  0/1:0,13 0/0:0,14 1/1:0,15
    1      977330  .  T   C   .    .      AF=0.05 GT:AD  0/1:0,16 0/0:0,17 ./.:0,18
    1      977516  .  -   C   .    .      AF=0.06 GT:AD  1/1:0,19 1/1:0,20 ./.:0,21
    1      977570  .  G   A   .    .      AF=0.07 GT:AD  1/1:0,22 1/1:0,23 ./.:0,24
    1      978604  .  CT  -   .    .      AF=0.08 GT:AD  1/1:0,25 1/1:0,26 ./.:0,27
    1      978628  .  C   T   .    .      AF=0.09 GT:AD  ./.:28,0 0/0:29,0 ./.:30,0

Open the VCF, making sure to carry over the `AF`, `GT`, and `AD` data::

    >>> vcf = VCFTable.from_vcf("test.vcf.gz", Genome("hg19"), info_ids=["AF"], fmt_ids=["GT", "AD"])
    >>> vcf
    <VCFTable, len() = 10>

    >>> vcf[0]
    <VCFVariant chr1:949522:C:T:hg19>

Get a variant's INFO attribute::

    >>> vcf[5].AF
    0.0500000007451

    >>> vcf.info("AF")
    array([ 0.  ,  0.01,  0.02,  ..., 0.07,  0.08,  0.09], dtype=float32)

Query an interval::

    >>> interval = vcf[5].expand(300, 300)   # chr1:977030-977630
    >>> variants = vcf.find_within(interval)
    >>> variants
    [<VCFVariant chr1:977329:T:C:hg19>,
     <VCFVariant chr1:977515::C:hg19>,
     <VCFVariant chr1:977569:G:A:hg19>]

Get indices of variants returned by a query::

    >>> indices = [vcf.index_of(v) for v in variants]
    >>> indices
    [5, 6, 7]

Get per-sample genotype and allelic depth for specific variants::

    >>> gt = vcf.format('GT')
    >>> gt.shape
    (10L, 3L)
    >>> gt[indices]
    array([[1, 0, 0],
           [2, 2, 0],
           [2, 2, 0]], dtype=int8)

    >>> ad = vcf.format('AD')
    >>> ad.shape
    (10L, 3L)
    >>> ad[indices]
    array([[[ 0, 16],
            [ 0, 17],
            [ 0, 18]],

           [[ 0, 19],
            [ 0, 20],
            [ 0, 21]],

           [[ 0, 22],
            [ 0, 23],
            [ 0, 24]]], dtype=int32)

The fastest way to filter variants by FORMAT columns is to use numpy on the entire array::

    >>> mask = np.any(ad.sum(axis=2) >= 20, axis=1)  # Find variants with at least one
    >>> variants = vcf.where(mask)                   # sample having ad >= 20
    >>> variants
    [<VCFVariant chr1:977515::C:hg19>,
     <VCFVariant chr1:977569:G:A:hg19>,
     <VCFVariant chr1:978603:CT::hg19>,
     <VCFVariant chr1:978627:C:T:hg19>]

You can run the above examples with ``$python demos/query_vcf.py`` from the GenomeKit directory.

Feature Extraction with Variants
--------------------------------

Besides reference genomes (:py:class:`.Genome`), GenomeKit can also represent
a genome that has variants applied to it (:py:class:`.VariantGenome`).
The idea is that you can pass either a reference genome or a variant genome
to your feature extraction code, and both cases will work transparently.

:py:class:`~.VariantGenome` will accept either a single :py:class:`.Variant` object,
a variant string, or a list which can contain a mixture of both. When a `Variant` object
is passed it is checked that it is defined on the same reference genome as the `VariantGenome`.
When a string is supplied, as `Variant` object is created on the same reference genome as
the `VariantGenome`.

Consider a toy example, where the only feature we extract is
the DNA sequence flanking the 5' end of a CFTR transcript::

    def extract_features(genome):
        tran = genome.transcripts["ENST00000426809.1"]   # CFTR transcript
        span = tran.end5.expand(2, 5)                    # 7nt span at 5' end
        return genome.dna(span)                          # extract DNA

    ref = Genome("gencode.v19")
    variants = [Variant.from_string("chr7:117120149:A:G", ref),     # rs397508328
                Variant.from_string("chr7:117120151:G:T", ref)]     # rs397508657
    var = VariantGenome(ref, variants)
    print(extract_features(ref))
    print(extract_features(var))

The above outputs reference DNA, and a variant DNA::

    CCATGCA
    CCGTTCA

However, unless one plans to only support SNVs and not INDELs,
it is important to specify how each query interval should be aligned with
respect to insertions/deletions. Read on and learn about
:py:class:`.VariantGenome` and about the concept of "anchors".


Variant genomes
^^^^^^^^^^^^^^^

Variant genomes are made by applying a zero or more variants to a reference
genome via the :py:class:`.VariantGenome` class.
If a list of variants is given, they are `all` applied as if they were a
single complex variant::

    ref = Genome("hg19")
    var1 = VariantGenome(ref, ref.variant("chr7:117120188:A:T"))    # rs397508673 (A>T)
    var2 = VariantGenome(ref, ref.variant("chr7:117120190:A:-"))    # rs397508710 (delA)
    var3 = VariantGenome(ref, [ref.variant(x) for x in ["chr7:117120188:A:T",
                               "chr7:117120190:A:-"]])  # both variants together

.. note::
   Variants are currently specified using string format ``chromosome:position:ref:alt`` where the position is DNA1,
   just like in ClinVar or the UCSC browser (technically VCF 4.0/4.1 standard). Clinvar also has a special format to
   denote indels that includes the    preceding base as padding, such as in ``'chr7:117,231,993:TCT:T'``. GenomeKit
   can handle this format but it is optional. Therefore ``'chr7:117,231,994:CT:'`` would be equivalent to the
   previous variant. GenomeKit also allows either ``-``, or ``.`` to denote an empty *ref* or *alt* field
   (*e.g.* ``'chr7:117,231,994:CT:.'``). When the padding base in supplied, GenomeKit trims it off internally.
   GenomeKit also allows comma separators in the variant position (*e.g.* ``'117,231,993'``), as well as both
   UCSC (``'chr7'``) and ENSEMBL-style (``'7'``) chromosome names (only for nuclear chromosomes 1--23, X, and Y).

Given a query interval, the variant genome's
:py:meth:`~genome_kit.VariantGenome.dna`
method returns the variant sequence, rather than the reference sequence::

    >>> interval = Interval("chr7", "+", 117120185, 117120192, ref)
    >>> ref.dna(interval)
    'CCAAACT'
    >>> var1.dna(interval)   # (A>T)
    'CCTAACT'
    >>> var2.dna(interval)   # (delA)
    'CCAACT'
    >>> var3.dna(interval)   # (A>T, delA)
    'CCTACT'

Notice that when a length-changing variant falls within the query interval,
the length of the result changes. This is only the default behaviour.
The next section explains how to control alignment for feature extraction
by 'anchoring' an interval.


Length-changing variants
^^^^^^^^^^^^^^^^^^^^^^^^

Variants that insert or delete positions (INDELs) effectively change the
coordinate system of the variant genome. If an interval is specified on the
reference genome, and there are variants falling within that interval, then
the manner in which it should be lifted to the variant genome is ambiguous.
This section explains how you can control the lifting behaviour to
suit your feature extraction needs.

The mechanism that GenomeKit provides is `anchored` intervals.
The idea is that the user indicates a position within the reference interval
that remains aligned when the interval is lifted over to the variant genome.

For example, you may want to anchor your interval to its 5' or 3' end::

    >>> interval = Interval("chr7", "+", 117120185, 117120192, ref)
    >>> anchored_5p = interval.with_anchor("5p")  # Anchored to its 5' end
    >>> anchored_3p = interval.with_anchor("3p")  # Anchored to its 3' end

    >>> ref = Genome("hg19")
    >>> var = VariantGenome(ref, ref.variant("chr7:117120190:A:-"))  # rs397508710 (delA)
    >>> ref.dna(interval)
    'CCAAACT'
    >>> var.dna(interval)     # (shrink 3' end)
    'CCAACT'
    >>> var.dna(anchored_5p)  # (fill 3' end)
    'CCAACTT'
    >>> var.dna(anchored_3p)  # (fill 5' end)
    'TCCAACT'

Besides the `"5p"` and `"3p"` modes, there are other ways to anchor
your intervals.
See :ref:`anchors` for a more in-depth explanation.

Motif finding
-------------

GenomeKit can find motifs in both reference and mutant sequences for you. This
is based on string matching, not PWMs.

The example shows how to search for motifs on both the reference and a variant
genome using the :py:meth:`genome_kit.Genome.find_motif` and
:py:meth:`genome_kit.VariantGenome.find_motif`::

    >>> genome = Genome('hg19')

    >>> # Short sequence from CFTR
    >>> interval = Interval('chr7', '+', 117231957, 117232030, genome)
    >>> genome.dna(interval)
    'TTGATATTTATATGTTTTTATATCTTAAAGCTGTGTCTGTAAACTGATGGCTAACAAAACTAGGATTTTGGTC'

    >>> motif = 'AACAA'
    >>> matches = genome.find_motif(interval, motif)
    >>> matches
    [Interval("chr7", "+", 117232009, 117232009, "hg19", 117232009)]

The returned interval is empty but can be expanded
for feature extraction::

    >>> matches[0].expand(5, 5)
    Interval("chr7", "+", 117232004, 117232014, "hg19", 117232009)

The returned interval always has its anchor set equal to its position. This means
that the interval will always stay aligned
to the same position on the reference genome when variants are applied.

The default is to return the empty interval upstream of the
first nucleotide of the motif matches, *i.e.* the motif is
immediately downstream.
But in many cases we would like to change that behaviour. For example, when we
match the acceptor core splice site motif ``AG``, we usually want the match
to be aligned to the 3' end of the ``AG`` motif. For
this reason, :py:meth:`~genome_kit.Genome.find_motif` supports the
``match_position`` argument. The default is ``match_position=0`` (equivalently
``match_position='5p'``), which aligns the match to the 5' end of the motif.
In the above case of searching for acceptor motifs, we can set
``match_position='3p'`` to return potential splice sites::

    >>> motif = 'AG'
    >>> matches = genome.find_motif(interval, motif, match_position='3p')
    >>> matches
    [Interval("chr7", "+", 117231987, 117231987, "hg19", 117231987),
     Interval("chr7", "+", 117232020, 117232020, "hg19", 117232020)]

We can verify that the ``AG`` is now immediately upstream of the returned
empty intervals::

    >>> [genome.dna(match.expand(2, 0)) for match in matches]
    ['AG', 'AG']

Alternatively, ``match_position`` can also be an integer: ``match_position=0``
is equivalent to ``match_position='5p'``, ``match_position=len(motif)`` is
equivalent to ``match_position='3p'`` and integers within that range match
positions within the motif.

By default :py:meth:`~genome_kit.Genome.find_motif` returns only *non
overlapping* matches, but this can be configured with the
``find_overlapping_matches`` parameter::

    >>> interval = Interval('chr7', '+', 117231957, 117232030, genome)
    >>> motif = 'TT'
    >>> genome.dna(interval)
    'TTGATATTTATATGTTTTTA'
    >>> genome.find_motif(interval, motif, find_overlapping_motifs=False)
    [Interval("chr7", "+", 117231957, 117231957, "hg19", 117231957),
     Interval("chr7", "+", 117231963, 117231963, "hg19", 117231963),
     Interval("chr7", "+", 117231971, 117231971, "hg19", 117231971),
     Interval("chr7", "+", 117231973, 117231973, "hg19", 117231973),
     Interval("chr7", "+", 117231981, 117231981, "hg19", 117231981),
     Interval("chr7", "+", 117232022, 117232022, "hg19", 117232022),
     Interval("chr7", "+", 117232024, 117232024, "hg19", 117232024)]]
    >>> genome.find_motif(interval, motif, find_overlapping_motifs=True)
    [Interval("chr7", "+", 117231957, 117231957, "hg19", 117231957),
     Interval("chr7", "+", 117231963, 117231963, "hg19", 117231963),
     Interval("chr7", "+", 117231964, 117231964, "hg19", 117231964),
     Interval("chr7", "+", 117231971, 117231971, "hg19", 117231971),
     Interval("chr7", "+", 117231972, 117231972, "hg19", 117231972),
     Interval("chr7", "+", 117231973, 117231973, "hg19", 117231973),
     Interval("chr7", "+", 117231974, 117231974, "hg19", 117231974),
     Interval("chr7", "+", 117231981, 117231981, "hg19", 117231981),
     Interval("chr7", "+", 117232022, 117232022, "hg19", 117232022),
     Interval("chr7", "+", 117232023, 117232023, "hg19", 117232023),
     Interval("chr7", "+", 117232024, 117232024, "hg19", 117232024)]

The :py:meth:`~genome_kit.VariantGenome.find_motif` works the same on
:py:class:`~genome_kit.VariantGenome` as on the reference genome::

    >>> variant = genome.variant("chr7:117231980::TTAGTT")  # Insertion
    >>> variant_genome = VariantGenome(genome, variant)
    >>> motif = 'AG'
    >>> matches = variant_genome.find_motif(interval, motif,
    ...                                     match_position='3p')
    >>> matches
    [Interval("chr7", "+", 117231979, 117231979, "hg19", 117231979, 4),
     Interval("chr7", "+", 117231987, 117231987, "hg19", 117231987),
     Interval("chr7", "+", 117232020, 117232020, "hg19", 117232020)]

The result are the previous two matches plus a third one where we inserted
the motif into the reference sequence. This is a special case since the
*position within the insertion* has no alignment to the reference genome.
In this case the interval has its ``anchor_offset`` attribute set to four
to indicate that the motif matches at the fourth position in the insertion.


.. _sharing-data:

Sharing new data files
----------------------

Data such as assemblies, annotations, tracks, etc is stored in custom-built binary files.
APIs for building these files are provided as part of the API.
GenomeKit first searches for data files in the directory specified by the environment variable
``GENOMEKIT_DATA_DIR``, defaulting to ``appdirs.user_data_dir("genome_kit")``.
If the file is not found, by default GenomeKit attempts to download it from a public read-only
Google Cloud Storage bucket.

You can use your own GCS bucket by setting the environment variable ``GENOMEKIT_GCS_BUCKET``,
allowing you to upload and share your data files.

If you wish to use a different mechanism to store remote data files, you can provide an
implementation of :py:class:`genome_kit.data_manager.DataManager` and register it:

.. code-block:: python

    class MyDataManager(DataManager):
            def __init__(self, data_dir: str):
                ...

            def get_file(self, filename: str) -> str:
                ...

            def upload_file(self, filepath: str, filename: str, metadata: Dict[str, str]=None):
                ...

    gk.gk_data.data_manager = MyDataManager()

    # Alternatively, you can install a plugin package that advertises a DataManager implementation
    # endpoint under the "genomekit.plugins.data_manager" group. GenomeKit will automatically
    # use the plugin's data manager.
    # (see https://setuptools.pypa.io/en/latest/userguide/entry_point.html#entry-points-for-plugins)
