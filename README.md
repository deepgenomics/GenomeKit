# GenomeKit

[![run-unit-tests](https://github.com/deepgenomics/GenomeKit/actions/workflows/run-tests.yaml/badge.svg)](https://github.com/deepgenomics/GenomeKit/actions/workflows/run-tests.yaml)

## What is GenomeKit?

GenomeKit is Deep Genomics' Python library for fast and easy access to
genomic resources such as sequence, data tracks, and annotations.
The goal is to let machine learning researchers build data sets
easily, and to be creative about how those data sets are designed.

GenomeKit is also designed to work with genome variants, giving users
a powerful way to extract features for different genotypes.

In the future it will include other useful resources like
secondary structure, conservation, and more.

## Useful Features

- Fast querying of DNA sequence and genomic tracks.
- Fast querying and structured access for annotations (GENCODE, RefSeq, ...).
- Fast querying and structured access for short read alignments (SAM/BAM files).
- Fast querying and filtering for variants (VCF files).
- Interval and Variant objects that are convenient, lightweight, and standardized.
- Scan the genome for motifs.

## Resource Available

- Reference DNA: hg19, hg19.p13.plusMT, hg38, hg38.p12, hg38.p13, hg38.p14, mm10.p6, mm39, rn6, macFas5, susScr11.
  *  To add more, see [Genomes](docs/genomes.rst)
- Annotations: GENCODE, RefSeq, Ensembl, APPRIS.
- Tracks: conservation, RNA structure, methylation, nucleosome positions, ...

## Installation

### Install via conda

The best way to install GenomeKit is via the pre-compiled conda packages
available from [anaconda.org/conda-forge/genomekit](https://anaconda.org/conda-forge/genomekit).

You can install GenomeKit with

    conda install genomekit

## Full Documentation (including developer instructions)

https://deepgenomics.github.io/GenomeKit

## Acknowledgements

We would like to express our gratitude to the following individuals who have contributed to the development of GenomeKit:

    Steve Chan ([@s22chan](https://github.com/s22chan))
    Andrew Delong ([@andrewdelong](https://github.com/andrewdelong))
    Hannes Bretschneider ([@hannes-brt](https://github.com/hannes-brt))
    Khalid Zuberi ([@kzuberi](https://github.com/kzuberi))
    Albert Cui ([@albscui](https://github.com/albscui))
    Aaron Uhr ([@jwuhr123](https://github.com/jwuhr123))
    Peiqi (Mark) Wang ([@tt6746690](https://github.com/tt6746690))
    Tasso Karkanis ([@TassoKarkanis](https://github.com/TassoKarkanis))
    Troy Simpson ([@mrtron](https://github.com/mrtron))

Special thanks to the original author, Andrew Delong, for laying the foundation of this project,
and to Steve Chan for his continued support and contributions.

We appreciate the efforts and contributions of all past contributors, whose work has been invaluable to the growth and improvement of GenomeKit.
