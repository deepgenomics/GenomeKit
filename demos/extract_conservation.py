# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome

# Print the mean conservation of every annotated exon
genome = Genome("gencode.v19")
for exon in genome.exons:
    cons = genome.phastcons_mammal_46way(exon)
    print(exon, "\t", cons.mean())
