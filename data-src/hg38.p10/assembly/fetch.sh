#!/usr/bin/env bash

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz
faToTwoBit GRCh38.p10.genome.fa.gz hg38.p10.2bit
twoBitInfo hg38.p10.2bit stdout | sort -k2rn > hg38.p10.chrom.sizes