#!/usr/bin/env bash

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M15/GRCm38.p5.genome.fa.gz
faToTwoBit GRCm38.p5.genome.fa.gz mm10.p5.2bit
twoBitInfo mm10.p5.2bit stdout | sort -k2rn > mm10.p5.chrom.sizes