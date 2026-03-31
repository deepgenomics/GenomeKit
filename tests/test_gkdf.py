import importlib.util
import tempfile
import unittest
from pathlib import Path

from genome_kit import Genome, Interval
from genome_kit.df import from_parquet, to_parquet

from . import MiniGenome

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

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_intervals(self):
        intervals = [
            Interval("chr1", "+", 2000, 3000, "hg19"),
            Interval("chr4", "-", 5000, 6000, "hg19"),
        ]
        df = pl.DataFrame({"intervals": [intervals]}, schema={"intervals": pl.Object})

        path = self.tmp_dir_path / "list_of_intervals.parquet"
        to_parquet(df, path)
        re_df = from_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_genomes(self):
        genomes = [Genome("hg38.p12"), Genome("gencode.v41")]
        df = pl.DataFrame({"genomes": [genomes]}, schema={"genomes": pl.Object})

        path = self.tmp_dir_path / "list_of_genomes.parquet"
        to_parquet(df, path)
        re_df = from_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_transcripts(self):
        g = Genome("gencode.v41")
        transcripts = list(g.transcripts)[:10]
        df = pl.DataFrame(
            {"transcripts": [transcripts]}, schema={"transcripts": pl.Object}
        )

        path = self.tmp_dir_path / "list_of_transcripts.parquet"
        to_parquet(df, path)
        re_df = from_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_multiple_types(self):
        g = Genome("gencode.v41")

        interval = Interval("chr5", "+", 2000, 3000, "hg19")
        transcript = g.genes[0].transcripts[0]
        gene = g.genes[0]
        exon = g.exons[0]

        df = pl.DataFrame(
            {
                "interval": [interval],
                "transcript": [transcript],
                "gene": [gene],
                "exon": [exon],
            }
        )

        path = self.tmp_dir_path / "multiple_types.parquet"
        to_parquet(df, path)
        re_df = from_parquet(path, lazy=False)
        self.assertEqual(re_df["interval"].item(), df["interval"].item())
        self.assertEqual(re_df["transcript"].item(), df["transcript"].item())
        self.assertEqual(re_df["gene"].item(), df["gene"].item())
        self.assertEqual(re_df["exon"].item(), df["exon"].item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_multiple_genomes(self):
        # test dataframe with multiple reference genomes in a single column
        g1 = Genome("gencode.v41")
        g2 = Genome("ucsc_refseq.2017-06-25")

        genes = [g1.genes[0], g2.genes[0]]
        df = pl.DataFrame({"genes": genes}, schema={"genes": pl.Object})

        path = self.tmp_dir_path / "multiple_genomes.parquet"
        to_parquet(df, path)
        re_df = from_parquet(path, lazy=False)
        self.assertEqual(re_df["genes"][0], df["genes"][0])
        self.assertEqual(re_df["genes"][1], df["genes"][1])


if __name__ == "__main__":
    unittest.main()
