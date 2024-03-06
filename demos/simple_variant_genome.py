# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome
from genome_kit import Variant
from genome_kit import VariantGenome


def extract_features(genome):
    tran = genome.transcripts["ENST00000426809.1"]  # CFTR transcript
    span = tran.end5.expand(2, 5)  # 7nt span at 5' end
    return genome.dna(span)  # extract DNA


ref = Genome("gencode.v19")
var = VariantGenome(
    ref,
    [
        Variant.from_string(x, ref) for x in [
            'chr7:117120149:A:G',  # rs397508328
            'chr7:117120151:G:T'  # rs397508657
        ]
    ])

print(extract_features(ref))  # CCATGCA
print(extract_features(var))  # CCGTTCA
