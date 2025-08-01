.. _genomes:

=======
Genomes
=======

GenomeKit now supports building and using assemblies and annotations. For
assemblies, the schema will follow the UCSC format, and for annotations, they
must be specified in GENCODE/Ensembl/NCBI GFF3 formats.

Examples
--------

Clone the GenomeKit git repo to see scripts under data-src/ for examples of how to build annotation data files.

    .. code-block:: Bash

        git clone https://github.com/deepgenomics/GenomeKit.git
        pushd GenomeKit

Scripts under ``data-src`` are used to obtain and generate the data files:

- ``data-src/<assembly>/assembly`` for the assembly, e.g ``data-src/hg19/assembly``
- ``data-src/<assembly>/<annotation-source>/<annotation>``, e.g ``data-src/hg19/GENCODE/v26lift37``
-

Assemblies
^^^^^^^^^^

#. Generate a hash file

#. .. code-block:: Bash

     echo $(python -c 'import genome_kit as gk; print(gk.Genome._refg_hash("hg19"))') > hg19.hash

(replace ``hg19`` with the desired assembly name)

#. Copy the ``2bit``, ``chrom.sizes``, and ``chromAlias.txt`` files from
   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/.
   and the ``hash`` file you generated into your

   .. code-block:: Bash

     python -c 'import os ; import appdirs ; print(os.environ.get("GENOMEKIT_DATA_DIR", appdirs.user_data_dir("genome_kit")))'

   directory.

   If you need to generate from a `fasta`:

   #. .. code-block:: Bash

        conda create -n ucsc-tools ucsc-fatotwobit ucsc-twobitinfo
   #. .. code-block:: Bash

        conda activate ucsc-tools

   #. follow the instructions at https://genome.ucsc.edu/goldenPath/help/twoBit.html
   #. optionally create an `chromAlias.txt` with any contig aliases required.

Annotations
^^^^^^^^^^^

#. .. code-block:: Bash

     python -c 'import genome_kit as gk; print(gk.GenomeAnnotation.build_gencode("MY_ANNO.gff3", "MY_ANNO", gk.Genome("MY_ASSEMBLY")))'

#. Copy the resulting files into the

   .. code-block:: Bash

     python -c 'import appdirs; print(appdirs.user_data_dir("genome_kit"))'

   directory.

   The `.dganno` file contains the compiled ``GFF3`` and the `.cfg` file
   contains metadata, such as ``refg=hg38``.


Appris / MANE
-------------

When adding an annotation, you can also generate Appris/MANE data files for it if public data is available.

Appris files are available on https://apprisws.bioinfo.cnio.es/pub/releases/.
As noted in the `GenomeKit source code <https://github.com/deepgenomics/GenomeKit/blob/41c1fbe5c011d04504b08eef4435c50771c01471/genome_kit/_build_appris.py#L107>`_,
we maintain a partial archive of Appris.

Gencode
^^^^^^^

For Gencode annotations, first find the matching Ensembl release. For example, for Gencode v47, the matching Ensembl
release is 113. So you'll need to find an Appris release that includes e113. The earliest release that includes e113 is
2024_10.v49 (e113v49).

If our archive doesn't already includes this release, you'll need to download the release and add it to the archive.

RefSeq
^^^^^^

You can similarly find the matching release for RefSeq annotations. For example, for RefSeq v110, the earliest Appris
release you can find rs110 is 2023_05.v48 (rs110v48).

MANE
^^^^

For MANE releases, search through versions on https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/.
Each MANE release includes a `README_versions.txt` the releated Ensembl and RefSeq releases.

For help on building Appris/MANE:

  .. code-block:: Bash

    python -m genome_kit build --help
