import importlib.util

from genome_kit import Genome, Interval
from genome_kit.df import from_parquet, to_parquet
from . import MiniGenome

import tempfile
import unittest
from pathlib import Path

HAS_POLARS = importlib.util.find_spec("polars") is not None
if HAS_POLARS:
    import polars as pl

class TestGkdfRoundTrip(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp_dir.cleanup)
        cls.tmp_dir_path = Path(cls.tmp_dir.name)

    @unittest.skip("MiniGenome and Genome type mismatch")
    def test_genome(self):
        # plain reference genome as well as gencode and refseq annotations
        genomes = ["hg38.p12", "gencode.v41", "ucsc_refseq.2017-06-25"]

        for genome_str in genomes:
            g = MiniGenome(genome_str)
            df = pl.DataFrame({"genome": [g]})

            path = self.tmp_dir_path / f"{genome_str}.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)

            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_interval(self):
        interval = Interval("chr5", "+", 2000, 3000, "hg19")
        df = pl.DataFrame({"interval": [interval]})

        path = self.tmp_dir_path / "interval.parquet"
        to_parquet(df, path)
        re_df = from_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_transcript(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            transcript = g.genes[0].transcripts[0]
            df = pl.DataFrame({"transcript": [transcript]})

            path = self.tmp_dir_path / f"{genome_str}_transcript.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())


    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_gene(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            gene = g.genes[0]
            df = pl.DataFrame({"gene": [gene]})

            path = self.tmp_dir_path / f"{genome_str}_gene.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_exon(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            exon = g.exons[0]
            df = pl.DataFrame({"exon": [exon]})

            path = self.tmp_dir_path / f"{genome_str}_exon.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())


    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_intron(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            intron = g.introns[0]
            df = pl.DataFrame({"intron": [intron]})

            path = self.tmp_dir_path / f"{genome_str}_intron.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())


    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_cds(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            cds = g.cdss[0]
            df = pl.DataFrame({"cds": [cds]})

            path = self.tmp_dir_path / f"{genome_str}_cds.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())


    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_utr3(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            utr3 = g.utr3s[0]
            df = pl.DataFrame({"utr3": [utr3]})

            path = self.tmp_dir_path / f"{genome_str}_utr3.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())


    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_utr5(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            utr5 = g.utr5s[0]
            df = pl.DataFrame({"utr5": [utr5]})

            path = self.tmp_dir_path / f"{genome_str}_utr5.parquet"
            to_parquet(df, path)
            re_df = from_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())


    

if __name__ == "__main__":
    unittest.main()