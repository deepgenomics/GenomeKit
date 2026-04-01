from __future__ import annotations

import functools
import json
import warnings
from collections.abc import Callable
from inspect import signature
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import polars as pl

import genome_kit as gk
from genome_kit._optional import require_polars

from .gk_structs import CURRENT_VERSION, CellType, ColumnInfo, GkDfVersion
from .registry import GK_TO_STRUCT, get_registry


def _map_batches_safe(fn: Callable) -> Callable:
    """Helper function to wrap a UDF and run safely with polars map_batches.

    Polars has a bug in map_batches that incorrectly forwards the return_dtype argument
    to the UDF. See https://github.com/pola-rs/polars/issues/24840.

    Args:
        fn: The user defined function to wrap.

    Returns:
        A wrapped version of the UDF that can be safely used with map_batches.
    """
    sig = signature(fn)

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        accepted = sig.parameters
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in accepted}
        return fn(*args, **filtered_kwargs)

    return wrapper


def _detect_gk_cols(
    lf: pl.LazyFrame, infer_schema_length: int = 100
) -> dict[str, ColumnInfo]:
    """Detect columns in the LazyFrame that contains GenomeKit objects.

    Args:
        lf: The LazyFrame to inspect.
        infer_schema_length: The number of rows to use for schema inference when
            detecting GenomeKit columns.

    Returns:
        A dictionary mapping column names to the ColumnInfo dataclass containing the
        GkDfType and CellType for the column.
    """
    pl = require_polars()

    lf_cols = lf.collect_schema().names()

    target_cols = {}

    # datatype inference done on first n=infer_schema_length rows. Follows inference
    # logic from Polars DataFrames when rows are provided.
    # see https://github.com/pola-rs/polars/blob/1cd236c60c01572c5ec6fdd252d8b20218d7b440/py-polars/src/polars/dataframe/frame.py#L248-L251
    head = lf.head(infer_schema_length).collect()

    for col in lf_cols:
        # remove nulls for type inference, list/scalar cols depend on first non-null value
        vals = head.get_column(col).drop_nulls()  # removes scalar nulls

        # column only contains null values in the first infer_schema_length rows
        if len(vals) == 0:
            warnings.warn(
                f"Column {col} contains only null values in the first {infer_schema_length} rows, "
                "unable to infer type for serialization. Please ensure this column "
                "contains non-null values for accurate serialization."
            )
            continue

        first = vals[0]
        if type(first) == list:
            cell_type = CellType.LIST
            # ensure all values are lists within a col
            if not all(type(v) == list for v in vals):
                raise ValueError(
                    f"Column {col} contains mixed data types. Please ensure all cells are the same type before serialization."
                )
            # cannot use Polars list expressions since lists of GenomeKit objects are stored as pl.Object
            col_types = {type(item) for v in vals for item in v if item is not None}

        else:
            cell_type = CellType.SCALAR
            # ensure all values are not lists within a col
            if not all(type(v) != list for v in vals):
                raise ValueError(
                    f"Column {col} contains mixed data types. Please ensure all cells are the same type before serialization."
                )
            col_types = set(vals.map_elements(type, return_dtype=pl.Object))

        if len(col_types) != 1:
            raise ValueError(
                f"Column {col} contains mixed data types. Please ensure all cells are the same type before serialization."
            )
        col_type = GK_TO_STRUCT.get(col_types.pop(), None)

        if col_type is None:
            # column is not a genomekit type, so no serialization needed
            continue

        target_cols[col] = ColumnInfo(cell_type=cell_type, gkdf_type=col_type)

    return target_cols


def _list_serializer(
    serializer: Callable[[pl.Series], pl.Series], return_dtype: Any
) -> Callable[[pl.Series], pl.Series]:
    """Helper function to convert a serializer to accept lists of GenomeKit objects.

    Args:
        serializer: A serializer function for a series of GenomeKit objects
        return_dtype: The return data type for the serialized series

    Returns:
        A serializer function for a series of lists of GenomeKit objects.
    """
    pl = require_polars()

    def _serialize_list(s: pl.Series) -> pl.Series:
        return pl.Series(
            name=s.name,
            values=[
                serializer(pl.Series(values=l)).to_list() if l is not None else None
                for l in s
            ],
            dtype=return_dtype,
        )

    return _serialize_list


def _init_gk_annotations(
    lf: pl.LazyFrame, target_cols: dict[str, dict]
) -> list[gk.Genome]:
    """Initialize GenomeKit annotations for all unique genomes in the LazyFrame.

    Prevents race conditions when opening dganno files during polars operations.
    Objects are returned in a list to keep weak references alive.

    Args:
        lf: The LazyFrame containing the serialized GenomeKit objects.
        target_cols: A dictionary mapping column names to their column information.
            target_cols is a dict from the ColumnInfo dataclass.

    Returns:
        A list of initialized gene tables for the unique genomes in the LazyFrame.
    """
    pl = require_polars()

    annotations = []

    # extract genome_str field from every column
    genomes_exprs = []
    genomes_list_exprs = []
    for c in target_cols.keys():
        if target_cols[c]["cell_type"] == CellType.SCALAR:
            genomes_exprs.append(pl.col(c).struct.field("genome_str"))
        else:
            genomes_list_exprs.append(pl.col(c).explode().struct.field("genome_str"))

    # expressions to extract genome_str must be run separately since exploded lists
    # may have more rows than the original dataframe
    plans = []

    if genomes_exprs:
        plans.append(
            lf.select(
                pl.concat_list(genomes_exprs)
                .explode()
                .drop_nulls()
                .unique()
                .alias("genome_str")
            )
        )

    if genomes_list_exprs:
        plans.append(
            lf.select(
                pl.concat_list(genomes_list_exprs)
                .explode()
                .drop_nulls()
                .unique()
                .alias("genome_str")
            )
        )

    genomes = pl.concat(plans).unique().collect()["genome_str"].to_list()

    # warms annotations for all unique annotation genomes in the file.
    # all annotations available for serialization are contained in dganno file
    for genome_str in genomes:
        genome = gk.Genome(genome_str)
        if genome.config == genome.reference_genome:
            # identifies reference genomes, instead of annotation genomes
            continue
        annotations.append(gk.Genome(genome_str).genes)

    return annotations


def _validate_gkdf_metadata(metadata: dict[str, str]) -> None:
    """Validate the parquet metadata for a gkdf parquet file.

    Args:
        metadata: The parquet metadata to validate.
    """

    try:
        version = GkDfVersion(metadata.get("gkdf_version"))
        if version not in GkDfVersion:
            raise ValueError(
                f"Unrecognized gkdf version in Parquet metadata, expected one of {[v.value for v in GkDfVersion]}"
            )
    except ValueError:
        raise ValueError(
            "Invalid or missing gkdf_version in Parquet metadata, unable to deserialize GenomeKit objects. "
        )
    if metadata.get("target_cols") is None:
        raise ValueError(
            "Missing target_cols in Parquet metadata, unable to deserialize GenomeKit objects."
        )


def _list_deserializer(
    deserializer: Callable[[pl.Series], pl.Series],
) -> Callable[[pl.Series], pl.Series]:
    """Helper function to convert a deserializer to accept lists of serialized GenomeKit objects.

    Args:
        deserializer: A deserializer function for a series of serialized GenomeKit objects

    Returns:
        A deserializer function for a series of lists of serialized GenomeKit objects.
    """
    pl = require_polars()

    def _deserialize_list(s: pl.Series) -> pl.Series:
        return pl.Series(
            name=s.name,
            values=[
                deserializer(pl.Series(values=l)).to_list() if l is not None else None
                for l in s
            ],
            dtype=pl.Object,
        )

    return _deserialize_list


def _deserialize_gk_cols(
    lf: pl.LazyFrame, target_cols: dict[str, dict]
) -> pl.LazyFrame:
    """Deserialize columns containing GenomeKit objects.

    Args:
        lf: The LazyFrame containing the serialized GenomeKit objects.
        target_cols: A dictionary mapping column names to their corresponding ColumnInfo.

    Returns:
        A LazyFrame with deserialized GenomeKit objects in the target columns.
    """
    pl = require_polars()
    registry = get_registry()

    def _build_deserialization_expr(col: str) -> pl.Expr:
        col_info = target_cols[col]  # dict representation of ColumnInfo
        gkdf_type = col_info["gkdf_type"]
        if col_info["cell_type"] == CellType.LIST:
            deserializer = _list_deserializer(
                registry[CURRENT_VERSION][gkdf_type].deserializer
            )
        else:
            deserializer = registry[CURRENT_VERSION][gkdf_type].deserializer

        return (
            pl.col(col)
            .map_batches(
                _map_batches_safe(deserializer),
                return_dtype=pl.Object,
            )
            .alias(col)
        )

    # with_columns_seq provides a 2x speedup here over with_columns
    return lf.with_columns_seq(_build_deserialization_expr(col) for col in target_cols)


# TODO: add union of pd.DataFrame
def write_parquet(
    df: pl.DataFrame | pl.LazyFrame, path: str | Path, infer_schema_length: int = 100
) -> None:
    """Serialize a DataFrame with GenomeKit objects to a Parquet file.

    Args:
        df: A Polars DataFrame or LazyFrame with columns containing GenomeKit objects.
        path: The file path to write the Parquet file to.
        infer_schema_length: The number of rows to use for schema inference when writing the Parquet file.
    """
    pl = require_polars()

    path = Path(path)
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    # mapping from column name to ColumnInfo dataclass
    target_cols = _detect_gk_cols(df, infer_schema_length=infer_schema_length)

    if not target_cols:
        warnings.warn(
            "No GenomeKit columns detected for serialization, writing DataFrame as is."
        )
        df.sink_parquet(path)
        return

    registry = get_registry()

    def _build_serialization_expr(col: str) -> pl.Expr:
        col_info = target_cols[col]  # ColumnInfo dataclass
        gkdf_type = col_info.gkdf_type
        if col_info.cell_type == CellType.LIST:
            return_dtype = pl.List(inner=registry[CURRENT_VERSION][gkdf_type].struct)
            serializer = _list_serializer(
                registry[CURRENT_VERSION][gkdf_type].serializer,
                return_dtype=return_dtype,
            )
        else:
            return_dtype = registry[CURRENT_VERSION][gkdf_type].struct
            serializer = registry[CURRENT_VERSION][gkdf_type].serializer

        return (
            pl.col(col)
            .map_batches(
                _map_batches_safe(serializer),
                return_dtype=return_dtype,
            )
            .alias(col)
        )

    df = df.with_columns(_build_serialization_expr(col) for col in target_cols)

    # convert ColumnInfo dataclass to a serializable format
    target_col_metadata = {col: target_cols[col].to_dict() for col in target_cols}

    metadata = {
        "gkdf_version": CURRENT_VERSION.value,
        "gk_version": gk.__version__,
        "target_cols": json.dumps(target_col_metadata),
    }

    df.sink_parquet(path, metadata=metadata)


def read_parquet(path: str | Path, lazy: bool = False) -> pl.DataFrame | pl.LazyFrame:
    """Deserialize a Parquet file containing GenomeKit objects into a Polars DataFrame or LazyFrame.

    Args:
        path: The file path to read the Parquet file from.
        lazy: If True, return a LazyFrame. Otherwise, return a DataFrame.

    Returns:
        A Polars DataFrame or LazyFrame with deserialized GenomeKit objects.
    """
    pl = require_polars()

    path = Path(path)
    metadata = pl.read_parquet_metadata(path)
    _validate_gkdf_metadata(metadata)
    target_cols = json.loads(metadata.get("target_cols"))

    lf = pl.scan_parquet(path)

    # collect unique genome strings in the file and initialize, prevents race conditions
    # on opening dganno files.
    # genomes returned in dummy variable to keep weak reference alive for deserialization
    _ = _init_gk_annotations(lf, target_cols)

    lf = _deserialize_gk_cols(lf, target_cols)

    return lf if lazy else lf.collect()
