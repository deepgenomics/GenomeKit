=========
GenomeKit
=========

------------------
What is GenomeKit?
------------------

GenomeKit is Deep Genomics' Python library for fast and easy access to
genomic resources such as sequence, data tracks, and annotations.
The goal is to let machine learning researchers build data sets
easily, and to be creative about how those data sets are designed.

GenomeKit is also designed to work with genome variants, giving users
a powerful way to extract features for different genotypes.

In the future it will include other useful resources like
secondary structure, conservation, and more.

Useful Features
===============

- Fast querying of DNA sequence and genomic tracks.
- Fast querying and structured access for annotations (GENCODE, RefSeq, ...).
- Fast querying and structured access for short read alignments (SAM/BAM files).
- Fast querying and filtering for variants (VCF files).
- Interval and Variant objects that are convenient, lightweight, and standardized.
- Scan the genome for motifs.

Resources Available
===================

* Reference genomes:

  - hg19 (unpatched GRCh37)
  - hg19.p13.plusMT (patch 13 with rCRS MT)
  - hg38 (unpatched GRCh38)
  - hg38.p12
  - hg38.p13
  - hg38.p14
  - mm10.p6 (GRCm38)
  - mm39 (unpatched GRCm39)
  - rn6 (unpatched RGSC 6.0)
  - macFas5 (unpatched Macaca_fascicularis_5.0)
  - susScr11 (unpatched SGSC Sscrofa11.1)

* Annotations:

  - GENCODE v19, v29, v41, vM19, vM30, vM31 (basic; comprehensive)
  - UCSC RefSeq (snapshots)
  - NCBI RefSeq v105.20190906, v109, v110, m39.v109, Macfas_5.0.v101
  - ENSEMBL Rnor_6.0.88, Macfas_5.0.95, Sscrofa11.1.98

Next steps
==========

Read the :ref:`quickstart` tutorial to get started with GenomeKit and study the
:ref:`api` for details. GenomeKit also uses a convention called *anchors* to
disambiguate the alignment of reference and variant sequences that is explained
in detail in :ref:`anchors`.

Contents:

.. toctree::
    :maxdepth: 2

    quickstart
    anchors
    api
    genomes
    develop

------------------
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
