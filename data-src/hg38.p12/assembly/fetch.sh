#!/usr/bin/env bash

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.2bit
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chromAlias.txt
# missing alias; see https://groups.google.com/a/soe.ucsc.edu/g/genome/c/oXgnoLwXn1g/m/zLV4Wgb2AgAJ
if ! grep -q "chrUn_KI270752v1" hg38.p12.chromAlias.txt; then
    echo "chrUn_KI270752v1	HSCHRUN_RANDOM_CTG29	KI270752.1	NT_187507.1" >> hg38.p12.chromAlias.txt
fi
