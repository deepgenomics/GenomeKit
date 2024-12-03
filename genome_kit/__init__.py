# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

from importlib import metadata

from . import gk_data
from .data_manager import DataManager, DefaultDataManager, GCSDataManager
from .genome import Genome, ApprisNotAvailableError
from .genome_annotation import (
    Cds,
    CdsTable,
    Exon,
    ExonTable,
    Gene,
    GeneTable,
    GenomeAnnotation,
    Intron,
    IntronTable,
    Transcript,
    TranscriptTable,
    Utr,
    UtrTable,
)
from .genome_dna import GenomeDNA
from .genome_track import GenomeTrack, GenomeTrackBuilder
from .interval import Interval
from .jralign import (
    JReadAlignments,
    JunctionReadAlignment,
    JunctionReadAlignments,
    JunctionReadAlignmentsTable,
)
from .jrdist import (
    JunctionReadCount,
    JunctionReadDistribution,
    JunctionReadDistributionTable,
    ReadDistributions,
)
from .ralign import (
    Alignment,
    AlignmentMatch,
    AlignmentMatchTable,
    AlignmentTable,
    Junction,
    JunctionTable,
    ReadAlignments,
)
from .variant import Variant, VariantTable
from .variant_genome import VariantGenome
from .vcf_table import VCFTable, VCFVariant
from . import serialize

#########################################################################



#########################################################################


__version__ = metadata.distribution("genomekit").version


__all__ = [
    "__version__",
    "Alignment",
    "AlignmentMatch",
    "AlignmentMatchTable",
    "AlignmentTable",
    "Cds",
    "CdsTable",
    "DataManager",
    "Exon",
    "ExonTable",
    "Gene",
    "GeneTable",
    "Genome",
    "GenomeAnnotation",
    "GenomeDNA",
    "GenomeTrack",
    "GenomeTrackBuilder",
    "gk_data",
    "Interval",
    "Intron",
    "IntronTable",
    "JReadAlignments",
    "Junction",
    "JunctionReadAlignment",
    "JunctionReadAlignments",
    "JunctionReadAlignmentsTable",
    "JunctionReadCount",
    "JunctionReadDistribution",
    "JunctionReadDistributionTable",
    "JunctionTable",
    "ReadAlignments",
    "ReadDistributions",
    "Transcript",
    "TranscriptTable",
    "Utr",
    "UtrTable",
    "Variant",
    "VariantGenome",
    "VariantTable",
    "VCFTable",
    "VCFVariant",
]

#########################################################################
