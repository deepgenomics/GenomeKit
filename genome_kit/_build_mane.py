import os
import pickle

import genome_kit as gk
import csv
import gzip

from genome_kit import gk_data
from genome_kit._gk_data_config import get_mane_filename, get_mane_version

_PRIORITY = {
    "MANE Select": 0,
    "MANE Plus Clinical": 1,
}

_SUPPORTED_ANNOS = {
    "ncbi_refseq.hg38.p14_RS_2024_08": "1.4",
    "gencode.v41": "1.0",
}


def _unversioned_id(id: str) -> str:
    return id.rsplit(".", 1)[0]


def build_mane(mane_version: str, genome: gk.Genome):
    """Generates the data file that GK needs for fast retrieval of MANE Select transcripts.

    Parameters
    ----------
    mane_version : :str
        MANE version, e.g "1.4". see https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/
    genome : :py:class:~genome_kit:Genome
        A genome object with the matching annotations for the specified MANE version

    Returns
    -------
    tuple
        [1] a list of indices for :py:class:~genome_kit.Transcript objects for fast access to all
            the transcripts that have a matching MANE Select transcript
        [2] a list where the indices are that of the gene objects and the values are the index of the
        corresponding MANE Select transcript in the genome.transcripts list.
    """
    mane_transcript_by_gene_idx = [None] * len(genome.genes)

    # files originally downloaded from https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/
    mane_src = gk.gk_data.get_file(f"MANE.GRCh38.v{mane_version}.summary.txt.gz")
    with gzip.open(filename=mane_src, mode="rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        mane_status_idx = header.index("MANE_status")
        transcript_id_col_idx = (
            header.index("RefSeq_nuc")
            if "ncbi_refseq" in genome.config
            else header.index("Ensembl_nuc")
        )
        transcript_indices = set()

        for line in reader:
            # use only MANE Select (ignore MANE Plus Clinical)
            if line[mane_status_idx] != "MANE Select":
                continue

            mane_transcript_id = _unversioned_id(line[transcript_id_col_idx])
            try:
                transcript = genome.transcripts[mane_transcript_id]
            except KeyError:
                if os.environ.get("GENOMEKIT_QUIET") != "1":
                    print(
                        f"transcript {mane_transcript_id} not found in {genome.config}"
                    )
                continue

            transcript_idx = genome.transcripts.index_of(transcript)
            transcript_indices.add(transcript_idx)
            gene_idx = genome.genes.index_of(transcript.gene)
            assert (
                mane_transcript_by_gene_idx[gene_idx] is None
                or mane_transcript_by_gene_idx[gene_idx] == transcript_idx
            ), f"more than one MANE Select transcript for gene {transcript.gene.id}"
            mane_transcript_by_gene_idx[gene_idx] = transcript_idx

    return sorted(transcript_indices), mane_transcript_by_gene_idx


def build_full_mane_files(upload: bool = False):
    for anno in _SUPPORTED_ANNOS:
        genome = gk.Genome(anno)
        res = build_mane(get_mane_version(anno), genome)
        output_filename = get_mane_filename(anno)
        output_filepath = os.path.join(gk.gk_data._config["DATA_DIR"], output_filename)
        with open(output_filepath, "wb") as fp:
            pickle.dump(res, fp, 2)
        print("built", output_filepath)

        if upload is True:
            print(" uploading ...", end="")
            gk_data.upload_file(output_filepath, output_filename)


TEST_DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "tests", "data"
)
assert os.path.isdir(TEST_DATA_DIR)


def build_test_mane_file(mane_version, annotation):
    """Builds test GenomeKit MANE binary files, using mini-annotations for testing."""

    data_dir = os.path.join(TEST_DATA_DIR, "mini1")
    genome = gk.Genome(annotation)
    output = build_mane(mane_version, genome)
    output_filepath = os.path.join(data_dir, get_mane_filename(annotation))

    with open(output_filepath, "wb") as fp:
        pickle.dump(output, fp, 2)

    print("Wrote", output_filepath)


def build_test_mane_files():
    os.environ["GENOMEKIT_QUIET"] = "1"
    for anno, mane_ver in _SUPPORTED_ANNOS.items():
        build_test_mane_file(mane_ver, anno + ".mini")
