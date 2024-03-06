# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome

genome = Genome("gencode.v19")

# Get sorted list of unique exon intervals in the genome
intervals = sorted(set([exon.interval for exon in genome.exons]))

# Extract DNA for each interval and print it
for interval in intervals:
    dna = genome.dna(interval)
    if len(dna) > 30:
        dna = dna[:30] + '...'
    print("{}\t{}".format(interval.as_ucsc(), dna))
