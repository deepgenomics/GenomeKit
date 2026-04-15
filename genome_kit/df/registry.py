from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    import polars as pl

import genome_kit as gk
from genome_kit._optional import require_polars

from .gk_structs import GkDfType, GkDfVersion, get_structs

# mapping from GenomeKit object types to the gkdf type strings
GK_TO_STRUCT: dict[type[gk.GenomeAnnotation], GkDfType] = {
    gk.Genome: GkDfType.GENOME,
    gk.Interval: GkDfType.INTERVAL,
    gk.Transcript: GkDfType.TRANSCRIPT,
    gk.Gene: GkDfType.GENE,
    gk.Exon: GkDfType.EXON,
    gk.Intron: GkDfType.INTRON,
    gk.Cds: GkDfType.CDS,
    gk.Utr: GkDfType.UTR,
}


# entry for the gkdf registry
@dataclass
class GKTypeEntry:
    struct: pl.Struct
    serializer: Callable[[pl.Series], pl.Series]
    deserializer: Callable[[pl.Series], pl.Series]


_GKDF_TYPE_FIELD = "gkdf_type"
_SCHEMA_VERSION_FIELD = "schema_version"

SUPPORTED_VERSIONS = {v for v in GkDfVersion.__members__.values()}


@lru_cache(maxsize=1)  # cache to avoid recreating registry in same session
def get_registry() -> dict[GkDfVersion, dict[GkDfType, GKTypeEntry]]:
    """Fetch the registry containing serialization and deserilization functions.

    Returns:
        Dictionary mapping GkDfType to their corresponding serializer and deserializer
        functions, for each supported GkDfVersion.
    """
    pl = require_polars()
    gkdf_structs = get_structs()

    def _serialize_genome(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Genome objects by genome name."""
        return pl.Series(
            name=s.name,
            values=[
                (
                    {
                        _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
                        # config gives annotation genome name if applicable
                        "genome_name": genome.config,
                    }
                    if genome is not None
                    else None
                )
                for genome in s
            ],
            dtype=gkdf_structs[GkDfType.GENOME],
        )

    def _deserialize_genome(s: pl.Series) -> pl.Series:
        """Deserialize a Series of GenomeStruct back into GenomeKit Genome objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Genome(struct["genome_name"]) if struct is not None else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_interval(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Interval objects."""
        return pl.Series(
            name=s.name,
            values=[
                {
                    "chromosome": interval.chromosome,
                    "strand": interval.strand,
                    "start": interval.start,
                    "end": interval.end,
                    # intervals related to reference genome only
                    "refg": interval.reference_genome,
                }
                if interval is not None
                else None
                for interval in s
            ],
            dtype=gkdf_structs[GkDfType.INTERVAL],
        )

    def _deserialize_interval(s: pl.Series) -> pl.Series:
        """Deserialize a Series of IntervalStruct back into GenomeKit Interval objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Interval(
                    chromosome=struct["chromosome"],
                    strand=struct["strand"],
                    start=struct["start"],
                    end=struct["end"],
                    reference_genome=struct["refg"],
                )
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_transcript(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Transcript objects."""
        return pl.Series(
            name=s.name,
            values=[
                {
                    _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
                    "transcript_table_index": transcript.annotation_genome.transcripts.index_of(
                        transcript
                    ),
                    "anno": transcript.annotation_genome.config,
                }
                if transcript is not None
                else None
                for transcript in s
            ],
            dtype=gkdf_structs[GkDfType.TRANSCRIPT],
        )

    def _deserialize_transcript(s: pl.Series) -> pl.Series:
        """Deserialize a Series of TranscriptStruct back into GenomeKit Transcript objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Genome(struct["anno"]).transcripts[struct["transcript_table_index"]]
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_gene(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Gene objects."""
        return pl.Series(
            name=s.name,
            values=[
                {
                    _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
                    "gene_table_index": gene.annotation_genome.genes.index_of(gene),
                    "anno": gene.annotation_genome.config,
                }
                if gene is not None
                else None
                for gene in s
            ],
            dtype=gkdf_structs[GkDfType.GENE],
        )

    def _deserialize_gene(s: pl.Series) -> pl.Series:
        """Deserialize a Series of GeneStruct back into GenomeKit Gene objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Genome(struct["anno"]).genes[struct["gene_table_index"]]
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_exon(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Exon objects."""
        return pl.Series(
            name=s.name,
            values=[
                {
                    _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
                    "exon_table_index": exon.annotation_genome.exons.index_of(exon),
                    "anno": exon.annotation_genome.config,
                }
                if exon is not None
                else None
                for exon in s
            ],
            dtype=gkdf_structs[GkDfType.EXON],
        )

    def _deserialize_exon(s: pl.Series) -> pl.Series:
        """Deserialize a Series of ExonStruct back into GenomeKit Exon objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Genome(struct["anno"]).exons[struct["exon_table_index"]]
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_intron(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Intron objects."""
        return pl.Series(
            name=s.name,
            values=[
                {
                    _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
                    "intron_table_index": intron.annotation_genome.introns.index_of(
                        intron
                    ),
                    "anno": intron.annotation_genome.config,
                }
                if intron is not None
                else None
                for intron in s
            ],
            dtype=gkdf_structs[GkDfType.INTRON],
        )

    def _deserialize_intron(s: pl.Series) -> pl.Series:
        """Deserialize a Series of IntronStruct back into GenomeKit Intron objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Genome(struct["anno"]).introns[struct["intron_table_index"]]
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_cds(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Cds objects."""
        return pl.Series(
            name=s.name,
            values=[
                {
                    _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
                    "cds_table_index": cds.annotation_genome.cdss.index_of(cds),
                    "anno": cds.annotation_genome.config,
                }
                if cds is not None
                else None
                for cds in s
            ],
            dtype=gkdf_structs[GkDfType.CDS],
        )

    def _deserialize_cds(s: pl.Series) -> pl.Series:
        """Deserialize a Series of CDSStruct back into GenomeKit Cds objects."""
        return pl.Series(
            name=s.name,
            values=[
                gk.Genome(struct["anno"]).cdss[struct["cds_table_index"]]
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    def _serialize_utr(s: pl.Series) -> pl.Series:
        """Serialize a Series of GenomeKit Utr objects."""
        values = []
        for utr in s:
            if utr is None:
                values.append(None)
                continue
            ser_dict = {
                _SCHEMA_VERSION_FIELD: GkDfVersion.V1.value,
            }
            genome = utr.annotation_genome
            try:
                ser_dict["utr_table_index"] = genome.utr5s.index_of(utr)
                ser_dict["utr_type"] = "5prime"
            except ValueError:
                ser_dict["utr_table_index"] = genome.utr3s.index_of(utr)
                ser_dict["utr_type"] = "3prime"

            ser_dict["anno"] = genome.config
            values.append(ser_dict)

        # pl.Series constructor doesn't strictly enforce nested struct types. Errors
        # with "utr_type" field may silently manifest as null values.
        # TODO: check polars implementation, related issue #18841
        return pl.Series(
            name=s.name,
            values=values,
            dtype=gkdf_structs[GkDfType.UTR],
        )

    def _deserialize_utr(s: pl.Series) -> pl.Series:
        """Deserialize a Series of UtrStruct back into GenomeKit Utr objects."""
        return pl.Series(
            name=s.name,
            values=[
                (
                    gk.Genome(struct["anno"]).utr5s[struct["utr_table_index"]]
                    if struct["utr_type"] == "5prime"
                    else gk.Genome(struct["anno"]).utr3s[struct["utr_table_index"]]
                )
                if struct is not None
                else None
                for struct in s
            ],
            dtype=pl.Object,
        )

    REGISTRY: dict[GkDfVersion, dict[GkDfType, GKTypeEntry]] = {
        GkDfVersion.V1: {
            GkDfType.GENOME: GKTypeEntry(
                struct=gkdf_structs[GkDfType.GENOME],
                serializer=_serialize_genome,
                deserializer=_deserialize_genome,
            ),
            GkDfType.INTERVAL: GKTypeEntry(
                struct=gkdf_structs[GkDfType.INTERVAL],
                serializer=_serialize_interval,
                deserializer=_deserialize_interval,
            ),
            GkDfType.TRANSCRIPT: GKTypeEntry(
                struct=gkdf_structs[GkDfType.TRANSCRIPT],
                serializer=_serialize_transcript,
                deserializer=_deserialize_transcript,
            ),
            GkDfType.GENE: GKTypeEntry(
                struct=gkdf_structs[GkDfType.GENE],
                serializer=_serialize_gene,
                deserializer=_deserialize_gene,
            ),
            GkDfType.EXON: GKTypeEntry(
                struct=gkdf_structs[GkDfType.EXON],
                serializer=_serialize_exon,
                deserializer=_deserialize_exon,
            ),
            GkDfType.INTRON: GKTypeEntry(
                struct=gkdf_structs[GkDfType.INTRON],
                serializer=_serialize_intron,
                deserializer=_deserialize_intron,
            ),
            GkDfType.CDS: GKTypeEntry(
                struct=gkdf_structs[GkDfType.CDS],
                serializer=_serialize_cds,
                deserializer=_deserialize_cds,
            ),
            GkDfType.UTR: GKTypeEntry(
                struct=gkdf_structs[GkDfType.UTR],
                serializer=_serialize_utr,
                deserializer=_deserialize_utr,
            ),
        }
    }

    return REGISTRY
