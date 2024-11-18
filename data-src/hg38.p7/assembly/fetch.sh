#!/usr/bin/env bash

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz
faToTwoBit GRCh38.p7.genome.fa.gz hg38.p7.2bit
twoBitInfo hg38.p7.2bit stdout | sort -k2rn > hg38.p7.chrom.sizes

# need chromAlias for ncbi, using hg38.p12.chromAlias.txt
wget -O hg38.p7.chromAlias.txt https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chromAlias.txt
# missing alias; see https://groups.google.com/a/soe.ucsc.edu/g/genome/c/oXgnoLwXn1g/m/zLV4Wgb2AgAJ
if ! grep -q "KI270752.1" hg38.p7.chromAlias.txt; then
    echo "chrUn_KI270752v1	HSCHRUN_RANDOM_CTG29	KI270752.1	NT_187507.1" >> hg38.p7.chromAlias.txt
fi
