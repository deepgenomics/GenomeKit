# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from genome_kit import Genome


def has_coding_seq(e):
    return e.cds is not None  # Has CDS?


def has_good_level(e):
    return e.tran.level <= 2  # Is level 1 or 2?


def not_first_exon(e):
    return e.index > 0  # Is not first exon?


genome = Genome("gencode.v19")  # 1196293 exons
exons = filter(has_coding_seq, genome.exons)  # 724078 remaining
exons = filter(has_good_level, exons)  # 605573 remaining
exons = filter(not_first_exon, exons)  # 558099 remaining
sites = set(exon.end5 for exon in exons)  # 198992 unique

print("Found %d unique acceptor sites:" % len(sites))
for site in sites:
    print(genome.dna(site.expand(5, 5)))
