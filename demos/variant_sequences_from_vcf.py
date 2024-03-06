# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome
from genome_kit import VCFTable
from pprint import pprint
from tempfile import gettempdir


genome = Genome("hg19")

def print_combos(vcf):
    interval = genome.interval("chr1", "+", 949520, 949700)
    variants = vcf.sequence_variations(interval)
    combos = vcf.variant_combinations(variants)
    variant_sequences = [genome.variant_dna(interval, x) for x in combos]
    print("{} combinations:".format(len(combos)))
    pprint(combos)
    pprint(variant_sequences)


VCF = gettempdir() + "/variant_sequences_from_vcf.vcf"
with open(VCF, "w") as vcf:
    vcf.write("""##fileformat=VCFv4.3
##reference=GRCh37
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
1	949523	.	C	T	.	.	AF=0.50	GT	0/1
1	949608	.	G	A	.	.	AF=1.00	GT	1/1
1	949695	.	C	CG	.	.	AF=0.50	GT	1/.
""")

# Open the VCF, preserving AF and GT format fields.
vcf = VCFTable.from_vcf(VCF, genome, info_ids=["AF"], fmt_ids=["GT"], cache=False)
print_combos(vcf)

# Open the VCF, this time with low confidence call.
vcf = VCFTable.from_vcf(VCF, genome, info_ids=["AF"], fmt_ids={"GT": VCFTable.GT_HETEROZYGOUS}, cache=False)
print_combos(vcf)
