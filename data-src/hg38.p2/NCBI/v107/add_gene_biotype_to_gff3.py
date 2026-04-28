#!/usr/bin/env python3
"""
Add gene_biotype attribute to older NCBI RefSeq GFF3 files that lack it (e.g. v107).

The biotype is inferred from:
  - pseudo=true on gene lines -> "pseudogene"
  - gbkey on transcript lines -> mRNA means "protein_coding",
    ncRNA uses ncrna_class (lncRNA, miRNA, etc.), precursor_RNA -> "pre_miRNA",
    misc_RNA -> "misc_RNA"
  - If a gene has no transcripts, falls back to "misc_RNA"

Two-pass approach:
  1. First pass: collect transcript gbkey/ncrna_class per parent gene ID.
  2. Second pass: write out the file, injecting gene_biotype on gene lines.

Usage:
    python add_gene_biotype_to_gff3.py input.gff3.gz output.gff3.gz
    python add_gene_biotype_to_gff3.py input.gff3 output.gff3
"""

import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_attributes(attr_str):
    """Parse GFF3 attribute string into a dict."""
    attrs = {}
    for item in attr_str.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def infer_biotype_from_transcript(gbkey, ncrna_class):
    """Infer a gene_biotype string from transcript-level attributes."""
    if gbkey == "mRNA":
        return "protein_coding"
    if gbkey == "ncRNA":
        if ncrna_class:
            # ncrna_class values like lncRNA, miRNA, etc. map directly
            # to biotypes that GenomeKit knows (via as_biotype_compatible)
            return ncrna_class
        return "ncRNA"
    if gbkey == "precursor_RNA":
        return "pre_miRNA"
    if gbkey == "rRNA":
        return "rRNA"
    if gbkey == "tRNA":
        return "tRNA"
    # misc_RNA, other
    return "misc_RNA"


def main():
    parser = argparse.ArgumentParser(
        description="Add gene_biotype attribute to NCBI RefSeq GFF3 files that lack it."
    )
    parser.add_argument("input", help="Input GFF3 file (optionally gzipped)")
    parser.add_argument("output", help="Output GFF3 file (gzipped if .gz)")
    args = parser.parse_args()

    # --- Pass 1: collect transcript info per gene ---
    # gene_id -> list of (gbkey, ncrna_class)
    tran_info_by_gene = defaultdict(list)
    # gene_id -> bool (has pseudo=true)
    gene_is_pseudo = {}
    # gene_id -> bool (already has gene_biotype)
    gene_has_biotype = {}
    # Track which gene IDs exist
    gene_ids = set()

    print(f"Pass 1: scanning {args.input} ...", file=sys.stderr)
    with open_maybe_gzip(args.input) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                continue

            attrs = parse_attributes(cols[8])
            gff_id = attrs.get("ID", "")

            # Detect gene lines
            if cols[2] == "gene" or cols[2] == "pseudogene" or gff_id.startswith("gene"):
                gene_ids.add(gff_id)
                gene_is_pseudo[gff_id] = (
                    attrs.get("pseudo", "") == "true" or cols[2] == "pseudogene"
                )
                gene_has_biotype[gff_id] = "gene_biotype" in attrs

            # Detect transcript lines (rna IDs with Parent pointing to a gene)
            elif gff_id.startswith("rna"):
                parent = attrs.get("Parent", "")
                if parent in gene_ids or parent.startswith("gene"):
                    gbkey = attrs.get("gbkey", "")
                    ncrna_class = attrs.get("ncrna_class", "")
                    tran_info_by_gene[parent].append((gbkey, ncrna_class))

    # --- Compute biotypes ---
    gene_biotype = {}
    for gid in gene_ids:
        if gene_has_biotype[gid]:
            continue  # already present, don't touch
        if gene_is_pseudo[gid]:
            gene_biotype[gid] = "pseudogene"
        elif tran_info_by_gene[gid]:
            # Use the first transcript to infer
            gbkey, ncrna_class = tran_info_by_gene[gid][0]
            gene_biotype[gid] = infer_biotype_from_transcript(gbkey, ncrna_class)
        else:
            gene_biotype[gid] = "misc_RNA"

    n_added = len(gene_biotype)
    n_already = sum(1 for v in gene_has_biotype.values() if v)
    print(
        f"  {len(gene_ids)} genes found, {n_already} already have gene_biotype, "
        f"{n_added} will be added.",
        file=sys.stderr,
    )

    # --- Pass 2: write output with gene_biotype injected ---
    print(f"Pass 2: writing {args.output} ...", file=sys.stderr)
    with open_maybe_gzip(args.input) as fin, open_maybe_gzip(args.output, "wt") as fout:
        for line in fin:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                fout.write(line + "\n")
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                fout.write(line + "\n")
                continue

            attrs = parse_attributes(cols[8])
            gff_id = attrs.get("ID", "")

            if gff_id in gene_biotype:
                # Append gene_biotype to the attributes column
                cols[8] = cols[8].rstrip(";") + f";gene_biotype={gene_biotype[gff_id]}"
                fout.write("\t".join(cols) + "\n")
            else:
                fout.write(line + "\n")

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

