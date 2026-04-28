import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from genome_kit import Genome, Interval, Variant
from genome_kit.df import read_parquet, write_parquet
from genome_kit.df.gk_structs import CURRENT_VERSION

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
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)

            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_interval(self):
        interval = Interval("chr5", "+", 2000, 3000, "hg19")
        df = pl.DataFrame({"interval": [interval]})

        path = self.tmp_dir_path / "interval.parquet"
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_transcript(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            transcript = g.genes[0].transcripts[0]
            df = pl.DataFrame({"transcript": [transcript]})

            path = self.tmp_dir_path / f"{genome_str}_transcript.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_gene(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            gene = g.genes[0]
            df = pl.DataFrame({"gene": [gene]})

            path = self.tmp_dir_path / f"{genome_str}_gene.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_exon(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            exon = g.exons[0]
            df = pl.DataFrame({"exon": [exon]})

            path = self.tmp_dir_path / f"{genome_str}_exon.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_intron(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            intron = g.introns[0]
            df = pl.DataFrame({"intron": [intron]})

            path = self.tmp_dir_path / f"{genome_str}_intron.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_cds(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            cds = g.cdss[0]
            df = pl.DataFrame({"cds": [cds]})

            path = self.tmp_dir_path / f"{genome_str}_cds.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_utr3(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            utr3 = g.utr3s[0]
            df = pl.DataFrame({"utr3": [utr3]})

            path = self.tmp_dir_path / f"{genome_str}_utr3.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_utr5(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            utr5 = g.utr5s[0]
            df = pl.DataFrame({"utr5": [utr5]})

            path = self.tmp_dir_path / f"{genome_str}_utr5.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_variant(self):
        genomes = ["gencode.v41", "ucsc_refseq.2017-06-25"]
        for genome_str in genomes:
            g = Genome(genome_str)
            variant = Variant("chr1", 10000005, "G", "T", g)
            df = pl.DataFrame({"variant": [variant]})

            path = self.tmp_dir_path / f"{genome_str}_variant.parquet"
            write_parquet(df, path)
            re_df = read_parquet(path, lazy=False)
            self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_intervals(self):
        intervals = [
            Interval("chr1", "+", 2000, 3000, "hg19"),
            Interval("chr4", "-", 5000, 6000, "hg19"),
        ]
        df = pl.DataFrame({"intervals": [intervals]}, schema={"intervals": pl.Object})

        path = self.tmp_dir_path / "list_of_intervals.parquet"
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_genomes(self):
        genomes = [Genome("hg38.p12"), Genome("gencode.v41")]
        df = pl.DataFrame({"genomes": [genomes]}, schema={"genomes": pl.Object})

        path = self.tmp_dir_path / "list_of_genomes.parquet"
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_transcripts(self):
        g = Genome("gencode.v41")
        transcripts = list(g.transcripts)[:10]
        df = pl.DataFrame(
            {"transcripts": [transcripts]}, schema={"transcripts": pl.Object}
        )

        path = self.tmp_dir_path / "list_of_transcripts.parquet"
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_variants(self):
        g = Genome("gencode.v41")
        variants = [Variant("chr1", 10000005, "G", "T", g) for _ in range(10)]
        df = pl.DataFrame(
            {"variants": [variants]}, schema={"variants": pl.Object}
        )

        path = self.tmp_dir_path / "list_of_variants.parquet"
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
        self.assertEqual(re_df.item(), df.item())

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_list_of_gk_with_null(self):
        g = Genome("gencode.v41")
        transcripts = list(g.transcripts)[:10]
        transcripts[:3] = [None] * 3
        df = pl.DataFrame(
            {"transcripts": [transcripts]}, schema={"transcripts": pl.Object}
        )

        path = self.tmp_dir_path / "list_of_transcripts_with_null.parquet"
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
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
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
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
        write_parquet(df, path)
        re_df = read_parquet(path, lazy=False)
        self.assertEqual(re_df["genes"][0], df["genes"][0])
        self.assertEqual(re_df["genes"][1], df["genes"][1])

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_mismatch_types(self):
        # test that error is raised when cols have different types
        g = Genome("gencode.v41")
        gene = g.genes[0]
        interval = Interval("chr5", "+", 2000, 3000, "hg19")

        df = pl.DataFrame({"mixed": [gene, interval]}, schema={"mixed": pl.Object})
        path = self.tmp_dir_path / "mismatch_types.parquet"
        with self.assertRaises(ValueError):
            write_parquet(df, path)

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_mismatch_list_types(self):
        # test that error is raised when cols have different types
        g = Genome("gencode.v41")
        gene = g.genes[0]

        df = pl.DataFrame({"mixed": [gene, [gene]]}, schema={"mixed": pl.Object})
        path = self.tmp_dir_path / "mismatch_list_types.parquet"
        with self.assertRaises(ValueError):
            write_parquet(df, path)

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_no_gkdf_version(self):
        # test that error raised when no gkdf version is found in metadata
        df = pl.DataFrame({"genome": ["hg38.p12"]})

        path = self.tmp_dir_path / "no_gkdf_version.parquet"
        df.write_parquet(path, metadata={"some_other_key": "value"})
        with self.assertRaises(ValueError):
            read_parquet(path, lazy=False)

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_no_target_cols(self):
        # test that error raised when no target_cols is found in metadata
        df = pl.DataFrame({"genome": ["hg38.p12"]})

        path = self.tmp_dir_path / "no_target_cols.parquet"
        df.write_parquet(path, metadata={"gkdf_version": CURRENT_VERSION})
        with self.assertRaises(ValueError):
            read_parquet(path, lazy=False)

    @unittest.skipUnless(HAS_POLARS, "Polars is required for this genome_kit.df tests")
    def test_no_gk_version(self):
        # test that error raised when no gk version is found in metadata
        df = pl.DataFrame({"genome": ["hg38.p12"]})

        path = self.tmp_dir_path / "no_gk_version.parquet"
        target_cols = {"genome": {"cell_type": "scalar", "gkdf_type": "genome"}}
        df.write_parquet(
            path,
            metadata={
                "gkdf_version": CURRENT_VERSION,
                "target_cols": json.dumps(target_cols),
            },
        )
        with self.assertRaises(ValueError):
            read_parquet(path, lazy=False)


if __name__ == "__main__":
    unittest.main()
