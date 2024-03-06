# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome

genome = Genome("gencode.v19")

for gene in genome.genes:
    print(gene)
    for tran in gene.transcripts:
        print("  ", tran)
        for exon in tran.exons:
            print("     ", exon)
            if exon.cds:
                print("      ", exon.cds)
            if exon.utr5:
                print("      ", exon.utr5)
            if exon.utr3:
                print("      ", exon.utr3)
