# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
#
# This file allows the genome_kit package to be run as a script.
# The script provides subcommands related to GenomeKit development
# and maintenance, such as building internal data files.
#
from __future__ import annotations

import gzip
import os
import re
import tempfile
from contextlib import suppress
from pathlib import Path

from . import (
    Genome,
    GenomeAnnotation,
    _twobitutils,
    _util,
    gk_data,
)

# Which gencode annotation to use for which test genome
GENCODE_OR_NCBI_TEST_ANNOTATIONS = {
    "hg19.p13.plusMT": "gencode.v29lift37",  # v29lift37 for hg19
    "hg38.p12": "gencode.v29",  # v29 for hg38
    "hg38.p13": "gencode.v41",
    "hg38.p14": "ncbi_refseq.hg38.p14_RS_2024_08",
}

UCSC_REFSEQ_TEST_ANNOTATIONS = {
    "hg19": "ucsc_refseq.2017-06-25",
}


def _gz_open(filename):
    """Calls open() or gzip.open() depending on file extension"""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "rt")


def _check_file(filename):
    """Checks that a file exists in _gk_data, or the non-gzipped version exists.
    Returns the path to the .gz file, or to the non-gz file if the
    .gz file doesn't exist. If neither file exists, raises an error.
    """
    return gk_data.get_file(filename)


def _strip_gz(path):
    if path.endswith(".gz"):
        return path[:-3]
    return path


def _ensure_dir(path):
    _util.makedirs(path)
    return path


def _ensure_path(path):
    _ensure_dir(os.path.dirname(path))
    return path


def _tests_data_dir():
    """Returns a path to the GenomeKit/tests module in
    the GenomeKit source tree."""

    # Get patch to the genome_kit source tree directory.
    gk_root = os.path.dirname(os.path.dirname(__file__))
    gk_tests = os.path.join(gk_root, "tests")

    # Check if directory was found. If not, it probably means we were
    # called from an installed version of genome_kit.
    if not os.path.isdir(gk_tests):
        raise RuntimeError("Cannot find GenomeKit/tests directory. " "Can only build test files from develop mode.")

    gk_tests_data = os.path.join(gk_tests, "data")
    return gk_tests_data


##############################################################


class _GencodeBuilder(object):
    def __init__(self, refg, gff3_rel_path):
        self.refg = refg
        self.srcurl = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{}/{}".format(
            "mouse" if gff3_rel_path.startswith("release_M") else "human", gff3_rel_path)

    def fetch(self):
        # Simply download the file from original public URL. Only re-download if newer.
        return gk_data.wget(self.srcurl, timestamping=True, progress=True)

    def build(self, srcpath, dstpath):
        return GenomeAnnotation.build_gencode(srcpath, dstpath, Genome(self.refg))


class _NCBIBuilder(object):
    def __init__(self, refg, gff3_url):
        self.refg = refg
        self.srcurl = gff3_url

    def fetch(self):
        # Simply download the file from original public URL. Only re-download if newer.
        return gk_data.wget(self.srcurl, timestamping=True, progress=True)

    def build(self, srcpath, dstpath):
        return GenomeAnnotation.build_ncbi_refseq(srcpath, dstpath, Genome(self.refg))


class _UCSCRefSeqBuilder(object):
    def __init__(self, refg, version_str):
        self.refg = refg
        self.version_str = version_str

    def fetch(self):
        srcroot = Path(tempfile.gettempdir())
        srcpath = srcroot / f"genomekit.{self.version_str}"
        srcpath.mkdir(exist_ok=True)  # e.g. "/tmp/genomekit.ucsc_refseq.2017-06-25"

        for x in ["refGene.txt.gz", "refLink.txt.gz"]:
            resolved = gk_data.get_file(f"{self.version_str}/{x}")
            dst = srcpath / x
            if not dst.exists():
                dst.symlink_to(resolved)
        """ for latest
        gk_data.wget(url="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz",
                     dst=os.path.join(srcpath, "refGene.txt.gz"), timestamping=True, progress=True)
        gk_data.wget(url="http://hgdownload.cse.ucsc.edu/goldenpath/hgFixed/database/refLink.txt.gz",
                     dst=os.path.join(srcpath, "refLink.txt.gz"), timestamping=True, progress=True)
        """
        return str(srcpath)

    def build(self, srcpath, dstpath):
        return GenomeAnnotation.build_ucsc_refseq(srcpath, dstpath, Genome(self.refg))


_ANNOTATION_BUILDERS = {
    "gencode.v29":
        _GencodeBuilder("hg38.p12", "release_29/gencode.v29.annotation.gff3.gz"),
    "gencode.v41":
        _GencodeBuilder("hg38.p13", "release_41/gencode.v41.annotation.gff3.gz"),
    "gencode.v29lift37":
        _GencodeBuilder("hg19.p13.plusMT", "release_29/GRCh37_mapping/gencode.v29lift37.annotation.gff3.gz"),
    "ucsc_refseq.2017-06-25":
        _UCSCRefSeqBuilder("hg19", "hg19.ucsc_refseq.2017-06-25"),
    "ncbi_refseq.hg38.p14_RS_2024_08":
        _NCBIBuilder("hg38.p14", "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"),

}


##############################################################

# Genomic regions from which to extract test data. Must be contiguous per chrom
#
# MINI GENOME 1: chr2:74682101-74692600 (hg19, 1-based)
#                chr2:74454974-74465473 (hg38, 1-based)
#   Genes INO80B (9 gencode transcripts, + strand)
#         INO80B-WBP1 (2 gencode transcripts; fusion of INO80B and WBP1)
#         WBP1 (14 gencode transcripts; + strand)
#         MOGS (9 gencode transcripts; - strand; 15 clinvar variants;
#               PDB domains & structure)
#         SULT1A3 (mapped to multiple regions in chr16 for UCSC)
#   and they all fall within a 10,500nt window.
#
# NOTE: when changing, python -m genome_kit --build-2bit --build-anno --build-appris
#       must be run
TEST_GENOME_REGIONS = (
    ('mini1', 'hg19', [('chr1', 84971983, 85022178), ('chr2', 74682100, 74692599), ('chr16', 29471206, 30215650)], {}),
    ('mini1', 'hg19.p13.plusMT', [('chr1', 84971983, 85022178), ('chr2', 74682100, 74692599), ('chr16', 29471206, 30215650)], {}),
    ('mini1', 'hg38.p12', [('chr2', 74454973, 74465472)], {}),
    ('mini1', 'hg38.p13', [('chr2', 74454973, 74465472)], {}),
    ('mini1', 'hg38.p14', [('chr2', 74455022, 74465382)], {"chr2": "NC_000002.12"}), # coordinates adjusted for ncbi_refseq.hg38.p14_RS_2024_08
)


def build_test_2bit_file(name, refg, regions, chrom_aliases):
    """Extracts chrom:start-end (in 0-based coordinates, inclusive)
    from a 2bit file, creating a new 2bit file with only a single
    chromosome (chrom) made entirely of the extracted DNA
    (of length end-start+1).

    For example, if you extract chr5:5000-5999 from hg19.2bit,
    then the new 2bit file will contain a single chromosome named
    'chr5' containing 1000nt of sequence.
    """
    try:
        from twobitreader import TwoBitFile
    except ImportError:
        raise ImportError("To build test data files, follow the dev setup instructions to set up a dedicated environment")


    # Source 2bit file to excerpt
    infile = {  # not testing patch contigs yet, use standard assembly for convenience
        "hg19": _check_file("hg19.2bit"),
        "hg19.p13.plusMT": _check_file("hg19.2bit"),
        "hg38.p12": _check_file("hg38.2bit"),
        "hg38.p13": _check_file("hg38.2bit"), # for testing MANE 1.0 with gencode.v41
        "hg38.p14": _check_file("hg38.2bit"), # for testing MANE 1.4 with ncbi_refseq.hg38.p14_RS_2024_08
    }[refg]

    # Destination 2bit file (small)
    outdir = _tests_data_dir()
    outfile = _ensure_path(os.path.join(outdir, name, f"{refg}.mini.2bit"))

    # Build test 2bit by excerpting full-sized 2bit
    tb = TwoBitFile(infile)
    seqs = {}
    for chrom, start, end in regions:
        seqs[chrom] = tb[chrom][start:end + 1]
    _twobitutils.write2bit(outfile, seqs)
    Path(outfile).with_suffix(".chrom.sizes").write_text("\n".join(f"{chrom_aliases.get(chrom, chrom)} {end - start + 1}" for chrom, start, end in regions))
    print("Wrote", outfile)


def build_test_2bit_file_from_seq(name, chrom_dna_dict):
    """Builds a 2bit file from a { chromosome : dna_sequence } dict.
    The output file is always <name>/hg19.2bit because the tests.MiniGenome
    construct expects it."""
    outdir = _tests_data_dir()
    outfile = _ensure_path(os.path.join(outdir, "mini1", f"{name}.2bit"))
    _twobitutils.write2bit(outfile, chrom_dna_dict)
    Path(outfile).with_suffix(".chrom.sizes").write_text("\n".join(f"{k} {len(v)}" for k, v in chrom_dna_dict.items()))
    print("Wrote", outfile)


def build_test_2bit_files():
    """Builds the test 2bit files."""

    # Extract test genomic regions.
    for name, refg, regions, chrom_aliases in TEST_GENOME_REGIONS:
        build_test_2bit_file(name, refg, regions, chrom_aliases)

    #                                                            0         1         2         3         4
    # Build hand-specified 2bit files, used by some tests.       01234567890123456789012345678901234567890123
    build_test_2bit_file_from_seq('test_genome', {
        'chr1': 'NNNNACGTacgtACGTacgtACGTacgtACGTacgtACGTacgt',
        'chr2': 'AACCTTTTACGTAAACCCGGGTTTACCGGAAATGGATTAA'
    })
    build_test_2bit_file_from_seq('test_genome_alt', {'chr1': 'AAAA', 'chr2': 'CCCC'})


def build_test_annotation_gff3(regions, chrom_aliases, infile, outfile):
    """Extracts GFF3 elements contained entirely within chrom:start-end
    (in 0-based coordinates, inclusive) from a GFF3 file, creating a
    new GFF3 file.

    The new elements will have their coordinates changed to be relative
    to 'start', i.e. if the original 1-based coordinate in the GFF3 file
    was 5201 and 0-based coordinate start=5000, then the new 1-based
    coordinate in the GFF3 output will be 200.

    If an element falls within the interval but its parent element
    does not, it is excluded.
    """
    parent_ids = []
    id_re = re.compile('ID=([^;]+)')  # Find the ID of an element
    pa_re = re.compile('Parent=([^;]+)')  # Find the Parent ID of an element
    with open(outfile, "w") as out:
        with _gz_open(infile) as file:
            for line in file:
                if line.startswith("#"):
                    # Header/comment lines
                    if line.startswith("#description:"):
                        # Alter the "description" line
                        out.write(
                            line.rstrip()
                            + f" FILTERED BY {','.join(f'{chrom}:{start + 1}-{end + 1}' for chrom, start, end in regions)}\n"
                        )
                    elif line.startswith("##sequence-region"):
                        # Alter the "sequence-region" line that corresponds to this chromosome; omit others
                        with suppress(StopIteration):
                            chrom, start, end = next(
                                (chrom, start, end)
                                for chrom, start, end in regions
                                if line.startswith(f"##sequence-region {chrom} ")
                            )
                            out.write("##sequence-region %s 1 %d\n" % (chrom, end - start + 1))
                    else:
                        # Other comment lines get copied verbatim
                        out.write(line)
                else:
                    for chrom, start, end in (
                        (chrom, start, end)
                        for chrom, start, end in regions
                        if line.startswith(f"{chrom_aliases.get(chrom, chrom)}\t")
                    ):
                        # Element on the right chromosome, at least
                        lchrom, lsource, ltype, lstart, lend, lextra = line.split("\t", 5)
                        lstart = int(lstart) - 1
                        lend = int(lend)
                        if lstart >= start and lend <= end:
                            # Only add elements that either have no parent or whose parent
                            # has already been added because it too falls completely within
                            # the filter interval
                            parent_matches = pa_re.findall(lextra)
                            if ltype == "gene" or (len(parent_matches) > 0 and parent_matches[0] in parent_ids):
                                out.write("\t".join(
                                    [lchrom, lsource, ltype,
                                     str(lstart - start + 1),
                                     str(lend - start), lextra]))
                                # All elements below gene/transcript can only have gene/transcript parents IDs
                                # so those are the only ones we add to the list.
                                if ltype in ("gene", "transcript", "lnc_RNA", "mRNA"):
                                    parent_ids.append(id_re.findall(lextra)[0])


def build_test_ucsc_refseq_database(regions, srcdir, dstdir):
    """Extracts GFF3 elements contained entirely within chrom:start-end
    (in 0-based coordinates, inclusive) from a GFF3 file, creating a
    new GFF3 file.

    The new elements will have their coordinates changed to be relative
    to 'start', i.e. if the original 1-based coordinate in the GFF3 file
    was 5201 and 0-based coordinate start=5000, then the new 1-based
    coordinate in the GFF3 output will be 200.

    If an element falls within the interval but its parent element
    does not, it is excluded.
    """
    # Build refGene filtered by chrom and start-end interval
    transcript_ids = set()
    with open(os.path.join(dstdir, "refGene.txt"), "w") as out:
        with _gz_open(os.path.join(srcdir, "refGene.txt.gz")) as file:
            for line in file:
                lbin, transcript_id, lchrom, lstrand, lstart, lend, cds_start, cds_end, \
                    exon_count, exon_starts, exon_ends, lextra = line.split("\t", 11)
                lstart = int(lstart)
                lend = int(lend)
                for start in (
                    start
                    for chrom, start, end in regions
                    if lchrom == chrom and lstart >= start and lend <= end
                ):
                    assert exon_starts.endswith(",") and exon_ends.endswith(
                        ",")  # Assumption about UCSC table formatting
                    out.write("\t".join([
                        lbin, transcript_id, lchrom, lstrand,
                        str(lstart - start),
                        str(lend - start),
                        str(int(cds_start) - start),
                        str(int(cds_end) - start), exon_count,
                        ",".join([str(int(es) - start) for es in exon_starts.split(",")[:-1]]) + ",",
                        ",".join([str(int(ee) - start) for ee in exon_ends.split(",")[:-1]]) + ",", lextra
                    ]))
                    transcript_ids.add(transcript_id)

    # Build refLink table filtered by rows that contain an observed transcript ID
    with open(os.path.join(dstdir, "refLink.txt"), "w") as out:
        with _gz_open(os.path.join(srcdir, "refLink.txt.gz")) as file:
            for line in file:
                _, _, transcript_id, _ = line.split("\t", 3)
                if transcript_id in transcript_ids:
                    out.write(line)


def build_test_gencode_or_ncbi_file(name, refg, regions, chrom_aliases):
    """Builds a test GENCODE annotation data file.
    """
    if refg not in GENCODE_OR_NCBI_TEST_ANNOTATIONS:
        return

    # Destination DGANNO file
    dstfile = GENCODE_OR_NCBI_TEST_ANNOTATIONS[refg]
    dstpath = _ensure_path(os.path.join(_tests_data_dir(), name, f"{dstfile}.mini"))

    # Destination GFF3 file (to be filtered)
    mini_srcfile = dstfile + ".annotation.gff3"
    mini_srcpath = _ensure_path(os.path.join(_tests_data_dir(), name, mini_srcfile))

    # Fetch the full-sized GFF3 file (unfiltered)
    builder = _ANNOTATION_BUILDERS[dstfile]
    builder.refg = f"{builder.refg}.mini"
    full_srcpath = builder.fetch()

    # Build the filtered GFF3
    build_test_annotation_gff3(regions, chrom_aliases, full_srcpath, mini_srcpath)
    print("Wrote", mini_srcpath)

    # Build the binary annotation file
    dstpaths = builder.build(mini_srcpath, dstpath)
    print("Wrote", dstpaths)


def build_test_ucsc_refseq_file(name, refg, regions):
    """Builds a test UCSC RefSeq annotation data file.
    """
    if refg not in UCSC_REFSEQ_TEST_ANNOTATIONS:
        return

    # Destination DGANNO file
    dstfile = UCSC_REFSEQ_TEST_ANNOTATIONS[refg]
    dstpath = _ensure_path(os.path.join(_tests_data_dir(), name, f"{dstfile}.mini"))

    # Destination database directory (to be filtered)
    mini_srcdir = _ensure_dir(os.path.join(_tests_data_dir(), name, dstfile))

    # Fetch the full-sized database files
    builder = _ANNOTATION_BUILDERS[dstfile]
    builder.refg = f"{builder.refg}.mini"
    full_srcdir = builder.fetch()

    # Build the filtered database
    build_test_ucsc_refseq_database(regions, full_srcdir, mini_srcdir)
    print("Wrote", os.path.join(mini_srcdir, "refGene.txt"))
    print("Wrote", os.path.join(mini_srcdir, "refLink.txt"))

    # Build the binary annotation file
    dstpath = builder.build(mini_srcdir, dstpath)
    print("Wrote", dstpath)


def build_test_annotation_files():
    """Builds all test annotation data files.
    """

    for name, refg, regions, chrom_aliases in TEST_GENOME_REGIONS:
        build_test_gencode_or_ncbi_file(name, refg, regions, chrom_aliases)
        build_test_ucsc_refseq_file(name, refg, regions)
