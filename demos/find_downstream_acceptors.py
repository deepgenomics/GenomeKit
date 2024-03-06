# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Interval
from genome_kit import Genome

genome = Genome("gencode.v19")

# We want to map each of these to an acceptor at most 100nt downstream
coords = [
    Interval.from_dna0_coord("chr1", "-", 91661, "hg19"),
    Interval.from_dna0_coord("chr1", "-", 169295, "hg19"),
    Interval.from_dna0_coord("chr1", "+", 320861, "hg19")
]

for coord in coords:
    window = coord.expand(0, 100)  # Find all exons with
    exons = genome.exons.find_5p_within(window)  # acceptor in window
    for exon in exons:
        print(coord.as_ucsc(), "-->", exon)
