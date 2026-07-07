#!/usr/bin/env python3
"""Promote UCSC chromosome names (chr1, chr2, ...) to canonical for an assembly.

GenomeKit treats the first column of chromAlias.txt as the canonical chromosome
name: it is what intervals print/serialize as and what genome.chromosomes lists.
Assemblies fetched from GCA hubs (e.g. macFas6) use the GenBank accession as the
first column, so chr1/chr2/... only exist as aliases. This rewrites the assembly's
chromAlias.txt and chrom.sizes in place so the UCSC names are canonical instead.

The previous canonical names are kept as a (now non-first) alias column, so the
.2bit (whose sequence names reference them) and existing user code still resolve.
The UCSC column is located by name in the header, so column order is irrelevant.
Rows with no UCSC name keep their original canonical name.

Usage: promote_ucsc_canonical.py <chromAlias.txt> <chrom.sizes>
"""
import sys

UCSC_SCHEME = "ucsc"


def promote(alias_path, sizes_path):
    with open(alias_path) as f:
        lines = f.read().splitlines()

    # header format: "# ucsc<TAB>genbank<TAB>..."
    # find the ucsc column
    schemes = [s.strip() for s in lines[0].lstrip("#").strip().split("\t")]
    if UCSC_SCHEME not in schemes:
        sys.exit(f"{alias_path}: no column named {UCSC_SCHEME!r} (have {schemes})")
    sc = schemes.index(UCSC_SCHEME)

    # New column order: the UCSC column first, then the rest in their original order.
    order = [sc] + [i for i in range(len(schemes)) if i != sc]

    rename = {}  # original canonical (first column) -> new canonical (UCSC name)
    out = ["# " + "\t".join(schemes[i] for i in order)]
    for line in lines[1:]:
        if not line:
            continue
        fields = line.split("\t")
        canonical = fields[sc] if fields[sc] else fields[0]
        rename[fields[0]] = canonical
        out.append("\t".join(fields[i] for i in order))

    with open(alias_path, "w") as f:
        f.write("\n".join(out) + "\n")

    # Rewrite chrom.sizes' first column using the same old -> new mapping.
    with open(sizes_path) as f:
        sizes = f.read().splitlines()
    with open(sizes_path, "w") as f:
        for line in sizes:
            if not line:
                continue
            name, _, size = line.partition("\t")
            f.write(f"{rename.get(name, name)}\t{size}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    promote(sys.argv[1], sys.argv[2])
