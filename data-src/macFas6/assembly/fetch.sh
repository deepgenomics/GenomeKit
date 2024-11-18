#!/usr/bin/env bash

# See https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_011100615.1/ for source

wget https://hgdownload.soe.ucsc.edu/hubs/GCA/011/100/615/GCA_011100615.1/GCA_011100615.1.2bit
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/011/100/615/GCA_011100615.1/GCA_011100615.1.chrom.sizes.txt
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/011/100/615/GCA_011100615.1/GCA_011100615.1.chromAlias.txt

mv GCA_011100615.1.2bit macFas6.2bit
mv GCA_011100615.1.chrom.sizes.txt macFas6.chrom.sizes
mv GCA_011100615.1.chromAlias.txt macFas6.chromAlias.txt
