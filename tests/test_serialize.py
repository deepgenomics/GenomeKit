
# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
import pickle
import os
import unittest

from genome_kit import Genome, Interval, Variant, VariantGenome
from tests import MiniGenome


class TestSerialize(unittest.TestCase):
    def test_serialize(self):
        # must be a Genome since MiniGenomes aren't pooled and we do address
        # comparisons or gene/tran/etc objects
        genome = Genome("ucsc_refseq.2017-06-25")
        objs = [
            Interval("chr1", "+", 10, 20, "hg19"),
            [
                Interval("chrX", "-", 10, 200, "hg19"),
                Interval("chr2", "-", 10, 20, "hg19", 20),
            ],
            Variant("chr1", 10, "A", "", "hg19"),
            genome,
        ]
        if 'CI' not in os.environ:
            # avoid downloading full dgannos on CI
            objs += [
                genome.genes[0],
                genome.trans[-1],
                genome.exons[len(genome.exons) // 2],
                genome.introns[1],
                genome.cdss[-2],
                genome.utr3s[-3],
                genome.utr5s[-4],
            ]

        anno_genome = MiniGenome("gencode.v29")
        objs.append(VariantGenome(anno_genome, Variant.from_string('chr2:20:A:T', anno_genome)))

        # Check that pickling / unpickling works
        stored = pickle.dumps(objs, pickle.HIGHEST_PROTOCOL)
        loaded = pickle.loads(stored)
        self.assertEqual(loaded, objs)


if __name__ == "__main__":
    unittest.main()
