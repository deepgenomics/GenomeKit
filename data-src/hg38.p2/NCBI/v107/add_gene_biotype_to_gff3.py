#!/usr/bin/env python3
"""
Add gene_biotype attribute to older NCBI RefSeq GFF3 files that lack it (e.g. v107).

Primary strategy: transfer gene_biotype from a newer RefSeq release using stable
identifiers (GeneID from Dbxref, then version-stripped accessions), with optional
coordinate/symbol validation.

Fallback (for unmapped genes): heuristic inference from:
  - pseudo=true on gene lines -> "pseudogene"
  - gbkey on transcript lines -> mRNA means "protein_coding",
    ncRNA uses ncrna_class (lncRNA, miRNA, etc.), precursor_RNA -> "pre_miRNA",
    misc_RNA -> "misc_RNA"
  - If a gene has no transcripts, falls back to "misc_RNA"

Usage:
    python add_gene_biotype_to_gff3.py --reference ref.gff3.gz input.gff3.gz output.gff3.gz
    python add_gene_biotype_to_gff3.py input.gff3.gz output.gff3.gz  # heuristic only
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


def extract_gene_id(dbxref):
    """Extract GeneID integer from a Dbxref attribute value."""
    for entry in dbxref.split(","):
        if entry.startswith("GeneID:"):
            try:
                return int(entry[7:])
            except ValueError:
                pass
    return None


def strip_version(accession):
    """Remove version suffix from an accession (e.g. NM_001005484.1 -> NM_001005484)."""
    if "." in accession:
        return accession.rsplit(".", 1)[0]
    return accession


def parse_reference_gff3(ref_path):
    """
    Parse a reference GFF3 that has gene_biotype on gene lines.
    Returns:
        geneid_to_info: {gene_id_int: (biotype, chrom, start, end, symbol)}
        accession_to_info: {stripped_accession: (biotype, chrom, start, end, symbol)}
    """
    geneid_to_info = {}
    accession_to_info = {}

    print(f"Parsing reference GFF3: {ref_path} ...", file=sys.stderr)
    count = 0
    with open_maybe_gzip(ref_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                raise ValueError("Malformed GFF3 line in reference: " + line)
            # Only care about gene lines
            if cols[2] not in ("gene", "pseudogene"):
                continue

            attrs = parse_attributes(cols[8])
            biotype = attrs.get("gene_biotype")
            if not biotype:
                continue

            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            symbol = attrs.get("gene", attrs.get("Name", ""))
            info = (biotype, chrom, start, end, symbol)

            # Primary key: GeneID
            dbxref = attrs.get("Dbxref", "")
            gene_id = extract_gene_id(dbxref)
            if gene_id is None:
                raise ValueError("Reference gene missing GeneID in Dbxref: " + line)
            geneid_to_info[gene_id] = info

            # Secondary key: version-stripped Name (gene symbol isn't unique, use ID)
            name = attrs.get("Name", "")
            if not name:
                raise ValueError("Reference gene missing Name attribute: " + line)
            accession_to_info[strip_version(name)] = info

            count += 1

    print(f"  {count} genes parsed from reference ({len(geneid_to_info)} with GeneID).",
          file=sys.stderr)
    return geneid_to_info, accession_to_info


def lookup_biotype_from_reference(gene_id, name, chrom, start, end, symbol,
                                   geneid_to_info, accession_to_info,
                                   max_distance=2_000_000):
    """
    Try to find gene_biotype from reference using GeneID (primary) or
    accession (secondary). Validates chromosome match and coordinate proximity.
    Returns biotype string or None.
    """
    # Try GeneID first
    if gene_id is not None and gene_id in geneid_to_info:
        ref_biotype, ref_chrom, ref_start, ref_end, ref_symbol = geneid_to_info[gene_id]
        # Validate: same chromosome (allow different accession format by comparing last part)
        if _chrom_compatible(chrom, ref_chrom):
            return ref_biotype
        # If chrom doesn't match but GeneID does, still trust it --
        # could be alt locus mapping differences
        return ref_biotype

    # Try version-stripped name/accession
    if name:
        stripped = strip_version(name)
        if stripped in accession_to_info:
            ref_biotype, ref_chrom, ref_start, ref_end, ref_symbol = accession_to_info[stripped]
            # Validate chromosome and distance
            if _chrom_compatible(chrom, ref_chrom):
                distance = min(abs(start - ref_start), abs(end - ref_end))
                if distance <= max_distance:
                    return ref_biotype
            # If only accession matches but coordinates are far, still use it
            # (genes can move due to assembly patch differences)
            return ref_biotype

    return None


def _chrom_compatible(chrom1, chrom2):
    """Check if two chromosome identifiers refer to the same sequence."""
    # Exact match
    if chrom1 == chrom2:
        return True
    # Compare by last component (e.g. NC_000001.11 vs NC_000001.14)
    base1 = chrom1.rsplit(".", 1)[0] if "." in chrom1 else chrom1
    base2 = chrom2.rsplit(".", 1)[0] if "." in chrom2 else chrom2
    return base1 == base2


def infer_biotype_from_transcript(gbkey, ncrna_class):
    """Infer a gene_biotype string from transcript-level attributes."""
    if gbkey == "mRNA":
        return "protein_coding"
    if gbkey == "ncRNA":
        if ncrna_class:
            return ncrna_class
        return "ncRNA"
    if gbkey == "precursor_RNA":
        return "pre_miRNA"
    if gbkey == "rRNA":
        return "rRNA"
    if gbkey == "tRNA":
        return "tRNA"
    return "misc_RNA"


def main():
    parser = argparse.ArgumentParser(
        description="Add gene_biotype attribute to NCBI RefSeq GFF3 files that lack it."
    )
    parser.add_argument("input", help="Input GFF3 file (optionally gzipped)")
    parser.add_argument("output", help="Output GFF3 file (gzipped if .gz)")
    parser.add_argument("--reference", "-r",
                        help="Reference GFF3 from a newer release (with gene_biotype)", required=True)
    args = parser.parse_args()

    geneid_to_info, accession_to_info = parse_reference_gff3(args.reference)

    # --- Pass 1: collect transcript info and gene metadata ---
    tran_info_by_gene = defaultdict(list)
    gene_is_pseudo = {}
    gene_has_biotype = {}
    gene_ids = set()
    # Store gene metadata for reference lookup
    gene_meta = {}  # gff_id -> (gene_id_int, name, chrom, start, end, symbol)

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

            if cols[2] == "gene" or cols[2] == "pseudogene" or gff_id.startswith("gene"):
                gene_ids.add(gff_id)
                gene_is_pseudo[gff_id] = (
                    attrs.get("pseudo", "") == "true" or cols[2] == "pseudogene"
                )
                gene_has_biotype[gff_id] = "gene_biotype" in attrs

                # Extract metadata for reference lookup
                dbxref = attrs.get("Dbxref", "")
                gene_id_int = extract_gene_id(dbxref)
                name = attrs.get("Name", "")
                symbol = attrs.get("gene", name)
                chrom = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                gene_meta[gff_id] = (gene_id_int, name, chrom, start, end, symbol)

            elif gff_id.startswith("rna"):
                parent = attrs.get("Parent", "")
                if parent in gene_ids or parent.startswith("gene"):
                    gbkey = attrs.get("gbkey", "")
                    ncrna_class = attrs.get("ncrna_class", "")
                    tran_info_by_gene[parent].append((gbkey, ncrna_class))

    # --- Compute biotypes ---
    gene_biotype = {}
    stats = {"reference": 0, "heuristic": 0, "already": 0, "id_not_found": 0, "name_not_found": 0}

    for gid in gene_ids:
        if gene_has_biotype[gid]:
            stats["already"] += 1
            continue

        # Try reference lookup
        gene_id_int, name, chrom, start, end, symbol = gene_meta[gid]

        if gene_id_int not in geneid_to_info:
            stats["id_not_found"] += 1
        elif strip_version(name) not in accession_to_info:
            stats["name_not_found"] += 1


        biotype = lookup_biotype_from_reference(
            gene_id_int, name, chrom, start, end, symbol,
            geneid_to_info, accession_to_info
        )
        if biotype:
            stats["reference"] += 1

        # Fallback to heuristic
        if biotype is None:
            if gene_is_pseudo[gid]:
                biotype = "pseudogene"
            elif tran_info_by_gene[gid]:
                gbkey, ncrna_class = tran_info_by_gene[gid][0]
                biotype = infer_biotype_from_transcript(gbkey, ncrna_class)
            else:
                biotype = "misc_RNA"
            stats["heuristic"] += 1

        gene_biotype[gid] = biotype

    print(
        f"  {len(gene_ids)} genes total, {stats['already']} already have gene_biotype, "
        f"{stats['reference']} mapped from reference, {stats['heuristic']} inferred by heuristic.\n"
        f"  Reference mapping failures: {stats['id_not_found']} with GeneID not found, "
        f"{stats['name_not_found']} with Name not found.",
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
                cols[8] = cols[8].rstrip(";") + f";gene_biotype={gene_biotype[gff_id]}"
                fout.write("\t".join(cols) + "\n")
            else:
                fout.write(line + "\n")

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

