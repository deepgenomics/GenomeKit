#!/usr/bin/env bash

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M11/GRCm38.p4.genome.fa.gz
faToTwoBit GRCm38.p4.genome.fa.gz mm10.p4.2bit
twoBitInfo mm10.p4.2bit stdout | sort -k2rn > mm10.p4.chrom.sizes
# no chromAlias.txt available for mm10.p4
wget -O mm10.p4.chromAlias.txt https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/p6/mm10.p6.chromAlias.txt
