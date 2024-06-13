#!/usr/bin/env bash

set -euox pipefail

TMP_DIR=$(mktemp -d)
pushd "${TMP_DIR}"
DATA_DIR=$(python -c 'import appdirs; import os; print(os.environ.get("GENOMEKIT_DATA_DIR", appdirs.user_data_dir("genome_kit")))')

function write_hash_file() {
  refg_name=$1
  hashval=$(python -c 'import genome_kit as gk; print(gk.Genome._refg_hash("${refg_name}"))')
  echo ${refg_name} > "${hashval}.hash"
}

# hg38.p12
# To use: genome_kit.Genome("hg38.p12")
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.2bit
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p12/hg38.p12.chromAlias.txt
write_hash_file "hg38.p12"
# chrUn_KI270752v1 missing from some hg38 chromAlias.txt data sources.
# It was dropped from the RefSeq release due to being derived likely from the human-hamster CHO cell line.
# See https://groups.google.com/a/soe.ucsc.edu/g/genome/c/oXgnoLwXn1g/m/zLV4Wgb2AgAJ for more details.
if ! grep -q "chrUn_KI270752v1" hg38.p12.chromAlias.txt; then
    echo "chrUn_KI270752v1	HSCHRUN_RANDOM_CTG29	KI270752.1	NT_187507.1" >> hg38.p12.chromAlias.txt
fi
mv ./* "${DATA_DIR}"

## Gencode v29
## To use: genome_kit.Genome("gencode.v29")
wget -O v29.gff3.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_gencode("v29.gff3.gz", "gencode.v29", gk.Genome("hg38.p12"))'
rm v29.gff3.gz

## NCBI RefSeq v109
## To use: genome_kit.Genome("ncbi_refseq.v109")
wget -O v109.gff3.gz http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/GFF/ref_GRCh38.p12_top_level.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_ncbi_refseq("v109.gff3.gz", "ncbi_refseq.v109", gk.Genome("hg38.p12"))'
rm v109.gff3.gz

mv ./* "${DATA_DIR}"

# Mus musculus (house mouse) mm39
# To use: genome_kit.Genome("mm39")
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromAlias.txt
write_hash_file "mm39"
mv ./* "${DATA_DIR}"

## Gencode vM31
## To use: genome_kit.Genome("gencode.VM31")
wget -O vM31.gff3.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_gencode("vM31.gff3.gz", "gencode.vM31", gk.Genome("mm39"))'
rm vM31.gff3.gz

## NCBI m39.v109
## To use: genome_kit.Genome("ncbi_refseq.m39.v109")
wget -O m39.v109.gff3.gz http://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/109/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_ncbi_refseq("m39.v109.gff3.gz", "ncbi_refseq.m39.v109", gk.Genome("mm39"))'
rm m39.v109.gff3.gz

mv ./* "${DATA_DIR}"

# Rattus norvegicus (Norway rat) rn6
# To use: genome_kit.Genome("rn6")
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.2bit
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.chrom.sizes
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.chromAlias.txt
write_hash_file "rn6"

## Ensembl Rnor_6.0.88
## To use: genome_kit.Genome("Rnor_6.0.88")
wget -O Rnor_6.0.88.gff3.gz http://ftp.ensembl.org/pub/release-88/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.88.gff3.gz
python -c 'import genome_kit as gk; gk.GenomeAnnotation.build_gencode("Rnor_6.0.88.gff3.gz", "Rnor_6.0.88", gk.Genome("rn6"))'

mv ./* "${DATA_DIR}"
