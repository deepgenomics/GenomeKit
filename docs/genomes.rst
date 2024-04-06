.. _genomes:

=======
Genomes
=======

GenomeKit now supports building and using assemblies and annotations. For
assemblies, the schema will follow the UCSC format, and for annotations, they
must be specified in GENCODE/Ensembl/NCBI GFF3 formats.

Examples
--------

``starter/build.sh`` is a bash script that builds a few assemblies and annotations to get you started.

Assemblies
^^^^^^^^^^

#. Copy the `2bit`, `chrom.sizes`, and `chromAlias.txt` files from
   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/. into your

   .. code-block:: Bash

     python -c 'import appdirs; print(appdirs.user_data_dir("genome_kit"))'

   directory.

   If you need to generate from a `fasta`:

   #. .. code-block:: Bash

        mamba create -n ucsc-tools ucsc-fatotwobit ucsc-twobitinfo
   #. .. code-block:: Bash

        mamba activate ucsc-tools

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

