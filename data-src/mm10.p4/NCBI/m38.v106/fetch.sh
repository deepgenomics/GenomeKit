#!/usr/bin/env bash

source ../../../functions.sh
wget_ncbi_mm10 m38.v106.tmp.gff3.gz https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Mus_musculus/ARCHIVE/ANNOTATION_RELEASE.106/GFF/ref_GRCm38.p4_top_level.gff3.gz

# NW_012132909.1, NW_004450262.1, NW_004058055.1 aren't in mm10.p6.chromAlias.txt, which we are using for mm10.p4.
# ok to remove it because it's an old assembly

# doing extra work of unzipping and gzipping again, and grepping multiple times, but this only takes a few seconds longer
zgrep -v NW_012132909.1 m38.v106.tmp.gff3.gz | grep -v NW_004450262.1 | grep -v NW_004058055.1 | gzip > m38.v106.gff3.gz
