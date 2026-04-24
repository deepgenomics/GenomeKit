#!/usr/bin/env bash

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/GRCh38.p2.genome.fa.gz
faToTwoBit GRCh38.p2.genome.fa.gz hg38.p2.2bit
twoBitInfo hg38.p2.2bit stdout | sort -k2rn > hg38.p2.chrom.sizes

# need chromAlias for ncbi, using hg38.p12.chromAlias.txt
wget -O hg38.p2.chromAlias.txt https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chromAlias.txt
