#!/usr/bin/env bash

wget -O v107.src.gff3.gz https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.107/GFF/ref_GRCh38.p2_top_level.gff3.gz

# Download a recent RefSeq release as reference for gene_biotype transfer
wget -O ref_latest.gff3.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz

# The v107 gff is missing gene_biotype; transfer from the newer release, with heuristic fallback
python add_gene_biotype_to_gff3.py --reference ref_latest.gff3.gz v107.src.gff3.gz v107.tmp.gff3.gz

# NT_187507.1 isn't in hg38.p12.chromAlias.txt, which we are using for hg38.p2.
# ok to remove it because NCBI says:
# > Record suppressed. This RefSeq record was removed because the sequence was determined to have originated from Chinese hamster.
# https://www.ncbi.nlm.nih.gov/nuccore/NT_187507.1?report=genbank

# doing extra work of unzipping and gzipping again but saves code duplication and only takes a few seconds longer
zgrep -v NT_187507.1 v107.tmp.gff3.gz | gzip > v107.gff3.gz
rm v107.tmp.gff3.gz
