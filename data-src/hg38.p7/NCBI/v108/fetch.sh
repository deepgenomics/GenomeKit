#!/usr/bin/env bash

source ../../../functions.sh
wget_ncbi_hg38 v108.tmp.gff3.gz https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.108/GFF/ref_GRCh38.p7_top_level.gff3.gz

# NT_187507.1 isn't in hg38.p12.chromAlias.txt, which we are using for hg38.p7.
# ok to remove it because NCBI says:
# > Record suppressed. This RefSeq record was removed because the sequence was determined to have originated from Chinese hamster.
# https://www.ncbi.nlm.nih.gov/nuccore/NT_187507.1?report=genbank

# doing extra work of unzipping and gzipping again but saves code duplication and only takes a few seconds longer
zgrep -v NT_187507.1 v108.tmp.gff3.gz | gzip > v108.gff3.gz
rm v108.tmp.gff3.gz
