#!/usr/bin/env bash

# See https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036323735.1/ for source.
# UCSC does not yet provide a goldenPath/rn8 mirror, so pull from the assembly hub.

wget https://hgdownload.soe.ucsc.edu/hubs/GCF/036/323/735/GCF_036323735.1/GCF_036323735.1.2bit
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/036/323/735/GCF_036323735.1/GCF_036323735.1.chrom.sizes.txt
wget https://hgdownload.soe.ucsc.edu/hubs/GCF/036/323/735/GCF_036323735.1/GCF_036323735.1.chromAlias.txt

mv GCF_036323735.1.2bit rn8.2bit
mv GCF_036323735.1.chrom.sizes.txt rn8.chrom.sizes
mv GCF_036323735.1.chromAlias.txt rn8.chromAlias.txt

# This GCF hub uses RefSeq accessions (e.g. NC_086019.1) as the canonical chromosome
# names. Promote the UCSC names (chr1, chr2, ...) to canonical so the assembly behaves
# like the others; the accessions are retained as an alias column (required for the 2bit file)
python3 ../../promote_ucsc_canonical.py rn8.chromAlias.txt rn8.chrom.sizes
