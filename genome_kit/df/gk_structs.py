from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING

from genome_kit._optional import require_polars

if TYPE_CHECKING:  # import polars for type checking
    import polars as pl


class GkDfType(StrEnum):
    GENOME = "genome"
    INTERVAL = "interval"
    TRANSCRIPT = "transcript"
    GENE = "gene"
    EXON = "exon"
    INTRON = "intron"
    CDS = "cds"
    UTR = "utr"


class GkDfVersion(StrEnum):
    V1 = "1.0"


CURRENT_VERSION = GkDfVersion.V1


def get_structs() -> dict[GkDfType, pl.Struct]:
    """Return a mapping of GkDfType to their corresponding Polars Struct definitions."""
    pl = require_polars()

    GenomeStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("genome_str", pl.Utf8),  # reference or annotation genome
        ]
    )

    IntervalStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("chromosome", pl.Utf8),
            pl.Field("strand", pl.Utf8),
            pl.Field("start", pl.Int32),
            pl.Field("end", pl.Int32),
            pl.Field("genome_str", pl.Utf8),  # reference or annotation genome
        ]
    )

    TranscriptStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            # index of transcript within annotation genome transcript table
            pl.Field("transcript_table_index", pl.Int32),
            pl.Field("genome_str", pl.Utf8),  # annotation genome
        ]
    )

    GeneStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("gene_table_index", pl.Int32),
            pl.Field("genome_str", pl.Utf8),  # annotation genome
        ]
    )

    ExonStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("exon_table_index", pl.Int32),
            pl.Field("genome_str", pl.Utf8),  # annotation genome
        ]
    )

    IntronStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("intron_table_index", pl.Int32),
            pl.Field("genome_str", pl.Utf8),  # annotation genome
        ]
    )

    CdsStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("cds_table_index", pl.Int32),
            pl.Field("genome_str", pl.Utf8),  # annotation genome
        ]
    )

    UtrStruct = pl.Struct(
        [
            pl.Field("gkdf_type", pl.Utf8),
            pl.Field("schema_version", pl.Utf8),
            pl.Field("utr_type", pl.Utf8),  # "5prime" or "3prime"
            pl.Field("utr_table_index", pl.Int64),
            pl.Field("genome_str", pl.Utf8),  # annotation genome
        ]
    )

    return {
        GkDfType.GENOME: GenomeStruct,
        GkDfType.INTERVAL: IntervalStruct,
        GkDfType.TRANSCRIPT: TranscriptStruct,
        GkDfType.GENE: GeneStruct,
        GkDfType.EXON: ExonStruct,
        GkDfType.INTRON: IntronStruct,
        GkDfType.CDS: CdsStruct,
        GkDfType.UTR: UtrStruct,
    }
