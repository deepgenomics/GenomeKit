import os
import pickle

import genome_kit as gk
import csv
import gzip

from genome_kit import gk_data
from genome_kit._gk_data_config import get_mane_filename

_PRIORITY = {
    "MANE Select": 0,
    "MANE Plus Clinical": 1,
}

_SUPPORTED_MANE_VERSIONS_BY_ANNO = {
    "ncbi_refseq.hg38.p14_RS_2024_08": "1.4",
    "ncbi_refseq.hg38.p14_RS_2024_08.mini": "1.4",
    "gencode.v41": "1.0",
    "gencode.v41.mini": "1.0",
    "gencode.v46": "1.3",
}

def get_mane_version(annotation):
    """
    Get the MANE version based on the annotation version.
    :raises ValueError: if there is no matching MANE version for the specified annotation
    """
    try:
        return _SUPPORTED_MANE_VERSIONS_BY_ANNO[annotation]
    except KeyError:
        raise ValueError(f"MANE not supported for annotation {annotation}")


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
        [0] a list of indices for :py:class:~genome_kit.Transcript objects for fast access to all
            the transcripts that have a matching MANE Select transcript
        [1] a list where the indices are that of the gene objects and the values are the index of the
        corresponding MANE Select transcript in the genome.transcripts list.
        [2] a list of indices for :py:class:~genome_kit.Transcript objects for fast access to all
            the transcripts that have a matching MANE Plus Clinical transcript
        [3] a list where the indices are that of the gene objects and the values are lists of the
            indices of the corresponding MANE Plus Clinical transcripts in the genome.transcripts list.
    """
    select_transcript_by_gene_idx = [None] * len(genome.genes)
    plus_clinical_transcript_by_gene_idx = [set() for _ in range(len(genome.genes))]

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
        select_transcript_indices = set()
        plus_clinical_transcript_indices = set()

        for line in reader:
            mane_status = line[mane_status_idx]
            assert mane_status in ["MANE Select", "MANE Plus Clinical"]

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
            gene_idx = genome.genes.index_of(transcript.gene)

            if mane_status == "MANE Select":
                select_transcript_indices.add(transcript_idx)
                assert (
                    select_transcript_by_gene_idx[gene_idx] is None
                    or select_transcript_by_gene_idx[gene_idx] == transcript_idx
                ), f"more than one MANE Select transcript for gene {transcript.gene.id}"
                select_transcript_by_gene_idx[gene_idx] = transcript_idx

            elif mane_status == "MANE Plus Clinical":
                plus_clinical_transcript_indices.add(transcript_idx)
                plus_clinical_transcript_by_gene_idx[gene_idx].add(transcript_idx)

    return (sorted(select_transcript_indices),
            select_transcript_by_gene_idx,
            sorted(plus_clinical_transcript_indices),
            [sorted(list(x)) for x in plus_clinical_transcript_by_gene_idx])


def build_full_mane_files(upload: bool = False):
    for anno in _SUPPORTED_MANE_VERSIONS_BY_ANNO:
        if anno.endswith(".mini"):
            continue
        genome = gk.Genome(anno)
        res = build_mane(get_mane_version(anno), genome)
        output_filename = get_mane_filename(anno)
        output_filepath = os.path.join(gk.gk_data._config["DATA_DIR"], output_filename)
        with open(output_filepath, "wb") as fp:
            pickle.dump(res, fp, 4)
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
    output = list(build_mane(mane_version, genome))
    # assign fake MANE Plus Clinical transcripts to some of the genes in the test ranges
    output[2] = output[0]
    plus_clinical_transcript_by_gene_idx = output[3]
    transcript_indices = output[0]
    for transcript_idx in transcript_indices:
        transcript = genome.transcripts[transcript_idx]
        gene = transcript.gene
        gene_idx = genome.genes.index_of(gene)
        plus_clinical_transcript_by_gene_idx[gene_idx].append(transcript_idx)
        if len(gene.transcripts) > 1:
            another_tx = next(tx for tx in gene.transcripts if tx != transcript)
            plus_clinical_transcript_by_gene_idx[gene_idx].append(genome.transcripts.index_of(another_tx))
            plus_clinical_transcript_by_gene_idx[gene_idx] = sorted(plus_clinical_transcript_by_gene_idx[gene_idx])
    output[3] = plus_clinical_transcript_by_gene_idx
    output_filepath = os.path.join(data_dir, get_mane_filename(annotation))

    with open(output_filepath, "wb") as fp:
        pickle.dump(output, fp, 4)

    print("Wrote", output_filepath)


def build_test_mane_files():
    os.environ["GENOMEKIT_QUIET"] = "1"
    for anno, mane_ver in _SUPPORTED_MANE_VERSIONS_BY_ANNO.items():
        build_test_mane_file(mane_ver, anno + ".mini")
