# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome
from genome_kit import VCFTable
from tempfile import gettempdir
import numpy as np

# Dump a demo VCF file
VCF = gettempdir() + "/query_vcf.vcf"
with open(VCF, "w") as vcf:
    vcf.write("""##fileformat=VCFv4.2
##reference=GRCh37
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1	sample2	sample3
1	949523	.	C	T	.	.	AF=0.00	GT:AD	0/0:0,1	0/1:0,2	0/0:0,3
1	949608	.	G	A	.	.	AF=0.01	GT:AD	0/0:0,4	0/1:0,5	0/0:0,6
1	949695	.	C	CG	.	.	AF=0.02	GT:AD	0/0:0,7	0/1:0,8	0/1:0,9
1	949739	.	G	TC	.	.	AF=0.03	GT:AD	0/1:0,10	0/0:0,11	1/1:0,12
1	977028	.	G	T	.	.	AF=0.04	GT:AD	0/1:0,13	0/0:0,14	1/1:0,15
1	977330	.	T	C	.	.	AF=0.05	GT:AD	0/1:0,16	0/0:0,17	./.:0,18
1	977515	.	A	C	.	.	AF=0.06	GT:AD	1/1:0,19	1/1:0,20	./.:0,21
1	977570	.	G	A	.	.	AF=0.07	GT:AD	1/1:0,22	1/1:0,23	./.:0,24
1	978603	.	ACT	A	.	.	AF=0.08	GT:AD	1/1:0,25	1/1:0,26	./.:0,27
1	978628	.	C	T	.	.	AF=0.09	GT:AD	./.:28,0	0/0:29,0	./.:30,0
""")

# Open the VCF and start querying
vcf = VCFTable.from_vcf(VCF, Genome("hg19"), info_ids=["AF"], fmt_ids=["GT", "AD"])

print(vcf)

# Get AF of a variant as float
print(vcf[5].AF)

# Get AF of all variants as numpy array
print(vcf.info("AF"))

# Get a window around one of the variants (chr1:977030-977630)
window = vcf[5].expand(300, 300)

# Three returned variants:
#   chr1:977030:T:C
#   chr1:977516::C
#   chr1:977570:G:A
variants = vcf.find_within(window)
print(variants)

# Index of each variant
indices = [vcf.index_of(v) for v in variants]

# Each variant's per-sample genotypes
#   [[1 0 0]
#    [2 2 0]
#    [2 2 0]]
gt = vcf.format('GT')
print(gt[indices])

# Each variant's per-sample allelic depth
# [[[ 0 16]
#   [ 0 17]
#   [ 0 18]]
#
#  [[ 0 19]
#   [ 0 20]
#   [ 0 21]]
#
#  [[ 0 22]
#   [ 0 23]
#   [ 0 24]]]
ad = vcf.format('AD')
print(ad[indices])

# Use numpy to quickly find variants with at least one sample
# above an allelic depth threshold.
mask = np.any(ad.sum(axis=2) >= 20, axis=1)
variants = vcf.where(mask)
print(variants)
