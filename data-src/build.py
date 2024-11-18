import argparse
import shutil
from pathlib import Path
import logging
import os
import subprocess

import genome_kit as gk

FETCH_FILENAME = "fetch.sh"
ASSEMBLY_DIRNAME = "assembly"
ENSEMBL_DIRNAME = "Ensembl"
GENCODE_DIRNAME = "GENCODE"
NCBI_DIRNAME = "NCBI"
UCSC_DIRNAME = "UCSC"
CWD = Path(__file__).parent

logging.basicConfig(level=logging.INFO)

ANNO_DIRNAME2BUILDER = {
    GENCODE_DIRNAME: gk.GenomeAnnotation.build_gencode,
    ENSEMBL_DIRNAME: gk.GenomeAnnotation.build_gencode,
    NCBI_DIRNAME: gk.GenomeAnnotation.build_ncbi_refseq
}

ANNO_DIRNAME2PREFIX = {
    GENCODE_DIRNAME: "gencode",
    ENSEMBL_DIRNAME: "ensembl",
    NCBI_DIRNAME: "ncbi_refseq"
}

def run(target_dir: Path):
    if ASSEMBLY_DIRNAME in args.build_target:
        # assembly
        assembly_path = Path(CWD,args.build_target, FETCH_FILENAME)
        assert len(assembly_path.parts) >= 3, "expected at least 3 parts in assembly file path"
        assembly_name = assembly_path.parts[-3]
        print(f"{assembly_name=}")

        subprocess.run(["bash", "-euxo", "pipefail", "fetch.sh"], cwd=Path(CWD, assembly_path.parent), check=True)

        twobit_filepath = assembly_path.with_name(f"{assembly_name}.2bit")
        print(f"moving {twobit_filepath} to {target_dir}")
        shutil.move(twobit_filepath, os.path.join(target_dir, twobit_filepath.name))

        chromsizes_filepath = assembly_path.with_name(f"{assembly_name}.chrom.sizes")
        print(f"moving {chromsizes_filepath} to {target_dir}")
        shutil.move(chromsizes_filepath, os.path.join(target_dir, chromsizes_filepath.name))

        chromalias_filepath = assembly_path.with_name(f"{assembly_name}.chromAlias.txt")
        if chromalias_filepath.exists():
            print(f"moving {chromalias_filepath} to {target_dir}")
            shutil.move(chromalias_filepath, os.path.join(target_dir, chromalias_filepath.name))

        refg_hash = gk.Genome._refg_hash(assembly_name)
        hash_filename = os.path.join(target_dir, f"{refg_hash}.hash")
        print(f"writing {hash_filename}")
        with open(hash_filename, "w") as f:
            f.write(f"{assembly_name}")

        return

    # annotation
    dirname = next(dirname for dirname in [GENCODE_DIRNAME, ENSEMBL_DIRNAME, NCBI_DIRNAME] if dirname in args.build_target)
    prefix, builder = ANNO_DIRNAME2PREFIX[dirname], ANNO_DIRNAME2BUILDER[dirname]

    annotation_path = Path(CWD, args.build_target, FETCH_FILENAME)
    print(f"{annotation_path=}")
    assert len(annotation_path.parts) >= 2, "expected at least 2 parts in annotation file path"
    annotation_name = annotation_path.parts[-2]

    subprocess.run(["bash", "-euxo", "pipefail", "fetch.sh"], cwd=Path(CWD, annotation_path.parent), check=True)
    assembly_name = args.build_target.split("/")[0]
    files = builder(
        str(annotation_path.with_name(f"{annotation_name}.gff3.gz")), f"{prefix}.{annotation_name}", gk.Genome(assembly_name)
    )
    for f in files:
        print(f"moving {f} to {target_dir}")
        shutil.move(f, os.path.join(target_dir, f))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build genomekit data")
    parser.add_argument("build_target", type=str, help="The assembly/annotation to build")
    parser.add_argument("target_dir", type=str, help="The target directory to build the data")
    args = parser.parse_args()
    run(Path(args.target_dir))
