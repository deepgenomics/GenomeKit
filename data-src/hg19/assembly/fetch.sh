#!/usr/bin/env bash

if [ ! -f hg19.2bit ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
fi
if [ ! -f hg19.chrom.sizes ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
fi
if [ ! -f hg19.chromAlias.txt ]; then
    # need to filter out MT, which only exists in p13.plusMT
    wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chromAlias.txt | grep -vF chrMT > hg19.chromAlias.txt
fi
