.. _data_org:

------------------------
Remote Data Organization
------------------------

In order to facilitate fast operations and queries over genomic data, GenomeKit
uses an internal optimized binary data format. GenomeKit provides APIs for converting
existing data formats such as gff3, vcf, sam, bam, bed, wig, etc into its own format.

In addition, various data and metadata support files are used. Out of the various formats,
only files related to genomic assemblies and annotations can be downloaded from freely
provided remote storage.
Other files (genomic tracks, indexed VCFs, etc) are maintained by users.

Available File Formats
======================

Formats not specific to GenomeKit:

- ``2bit`` (raw assembly data)
- ``chrom.sizes`` / ``chromAlias.txt`` (text format)

GenomeKit-specific data files:

- ``dganno``: genomic annotation data
- ``cfg``: annotation metadata, such as associated assembly name (text format)
- ``appris.*.pkl``: Python-pickled appris metadata
- ``mane.*.pkl``: Python-pickled MANE metadata

Data Access
===========

Assembly and annotation files are accessed transparently by GenomeKit when users
create Genome objects. For example, when users call
``gk.Genome("gencode.v29").dna(gk.Genome("gencode.v29").exons[0].interval)``,
GenomeKit downloads (if not already locally available) the files:

- ``gencode.v29.cfg`` (annotation configuration)
- ``hg38.p12.2bit`` (associated assembly)
- ``hg38.p12.chrom.sizes`` (list of chromosomes/contigs and aliases)
- ``gencode.v29.v22.dganno`` (indexed annotation features ; "v22" specifies the format version)

In the source bucket, files are organized in a flat structure and are easily
found by their name. The binary version is matched via library code.

Files are also directly accessible by using GenomeKit's data_manager API.
This API also allows users to configure their own buckets (to also allow uploads)
or a different implementation.

Downloads are automatically cached locally and so are only downloaded once.
They are lazily downloaded, so that only the required files are retrieved.
