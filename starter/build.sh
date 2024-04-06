#!/usr/bin/env bash

set -euox pipefail

TMP_DIR=$(mktemp -d)
pushd "${TMP_DIR}"
DATA_DIR=$(python -c 'import appdirs; import os; print(os.environ.get("GENOMEKIT_DATA_DIR", appdirs.user_data_dir("genome_kit")))')

# hg38.p12
# To use: genome_kit.Genome("hg38.p12")
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.2bit
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chromAlias.txt
# chrUn_KI270752v1 missing from some hg38 chromAlias.txt data sources
if ! grep -q "chrUn_KI270752v1" hg38.p12.chromAlias.txt; then
    echo "chrUn_KI270752v1	HSCHRUN_RANDOM_CTG29	KI270752.1	NT_187507.1" >> hg38.p12.chromAlias.txt
fi
mv ./* "${DATA_DIR}"

## Gencode v29
## To use: genome_kit.Genome("gencode.v29")
wget -O v29.gff3.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_gencode("v29.gff3.gz", "gencode.v29", gk.Genome("hg38.p12"))'
rm v29.gff3.gz

## NCBI RefSeq v109
## To use: genome_kit.Genome("ncbi_refseq.v109")
wget -O v109.gff3.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/GFF/ref_GRCh38.p12_top_level.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_ncbi_refseq("v109.gff3.gz", "ncbi_refseq.v109", gk.Genome("hg38.p12"))'
rm v109.gff3.gz

mv ./* "${DATA_DIR}"

# macFas5
# To use: genome_kit.Genome("macFas5")
wget https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/bigZips/macFas5.2bit
wget https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/bigZips/macFas5.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/macFas5/bigZips/macFas5.chromAlias.txt
mv ./* "${DATA_DIR}"

## Ensembl 0.95
## To use: genome_kit.Genome("ensembl.Macfas_5.0.95")
wget -O Macfas_5.0.95.gff3.gz ftp://ftp.ensembl.org/pub/release-95/gff3/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.95.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_gencode("Macfas_5.0.95.gff3.gz", "ensembl.Macfas_5.0.95", gk.Genome("macFas5"))'
rm Macfas_5.0.95.gff3.gz

## NCBI RefSeq v101
## To use: genome_kit.Genome("ncbi_refseq.Macfas_5.0.v101")
wget -O Macfas_5.0.v101.gff3.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_fascicularis/annotation_releases/101/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.gff.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_ncbi_refseq("Macfas_5.0.v101.gff3.gz", "ncbi_refseq.Macfas_5.0.v101", gk.Genome("macFas5"))'
rm Macfas_5.0.v101.gff3.gz

mv ./* "${DATA_DIR}"
