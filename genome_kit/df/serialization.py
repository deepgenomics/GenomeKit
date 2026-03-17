import functools
import json
import time
import warnings
from collections.abc import Callable
from inspect import signature

import polars as pl

import genome_kit as gk

from .gk_structs import CURRENT_VERSION, GkDfType, GkDfVersion
from .registry import GK_TO_STRUCT, REGISTRY


def _map_batches_safe(fn: Callable):
    """Helper function to wrap a UDF and run safely with polars map_batches.

    Polars has a bug in map_batches that incorrectly forwards the return_dtype argument
    to the UDF. See https://github.com/pola-rs/polars/issues/24840.

    Args:
        fn: The user defined function to wrap.
    """
    sig = signature(fn)

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        accepted = sig.parameters
        filtered_kwargs = {k: v for k, v in kwargs.items() if k in accepted}
        return fn(*args, **filtered_kwargs)

    return wrapper


def detect_gk_cols(
    lf: pl.LazyFrame, columns: list[str] | None = None
) -> dict[str, GkDfType]:
    """Detect columns in the LazyFrame that contains GenomeKit objects.

    Args:
        lf: The LazyFrame to inspect.
        columns: Optional list of column names to check. If None, all columns will be checked.

    Returns:
        A dictionary mapping column names to their corresponding GenomeKit types.
    """

    lf_cols = lf.collect_schema().names()

    if not columns:
        columns = lf_cols

    target_cols = {}

    # polars Struct inferred from first row, same behaviour as Polars
    # see https://docs.pola.rs/user-guide/expressions/structs/#inferring-the-data-type-struct-from-dictionaries
    # materialize the first row to check data types, need the exact type not pl.Object
    first_row = lf.head(1).collect()[0]

    for col in columns:
        if col not in lf_cols:
            raise ValueError(
                f"Column '{col}' not found in the DataFrame, please check the column names and try again."
            )
        # item from first row of the column
        col_type = GK_TO_STRUCT.get(type(first_row[col][0]), None)

        if col_type is None:
            # column is not a genomekit type, so no serialization needed
            pass
        else:
            target_cols[col] = col_type

    return target_cols


# TODO: add union of pd.DataFrame
def to_parquet(
    df: pl.DataFrame | pl.LazyFrame,
    path: str,
    columns: list[str] | None = None,
) -> None:
    """Serialize a DataFrame with GenomeKit objects to a Parquet file.

    Args:
        df: A Polars DataFrame or LazyFrame with columns containing GenomeKit objects.
        path: The file path to write the Parquet file to.
        columns: Optional list of column names to serialize. If None, all GenomeKit
            columns will be serialized.

    """
    if isinstance(df, pl.DataFrame):
        df = df.lazy()

    target_cols = detect_gk_cols(df, columns)

    if not target_cols:
        warnings.warn(
            "No GenomeKit columns detected for serialization, writing DataFrame as is."
        )
        df.sink_parquet(path)
        return

    df = df.with_columns(
        pl.col(col)
        .map_batches(
            _map_batches_safe(REGISTRY[CURRENT_VERSION][target_cols[col]].serializer),
            return_dtype=REGISTRY[CURRENT_VERSION][target_cols[col]].struct,
        )
        .alias(col)
        for col in target_cols
    )

    metadata = {
        "gkdf_version": CURRENT_VERSION.value,
        "gk_version": gk.__version__,
        "target_cols": json.dumps(target_cols),
    }

    df.sink_parquet(path, metadata=metadata)


def _init_gk_annotations(lf: pl.LazyFrame, target_cols: dict[str, GkDfType]) -> None:
    """Initialize GenomeKit annotations for all unique genomes in the LazyFrame.

    Prevents race conditions when opening dganno files during polars operations.

    Args:
        lf: The LazyFrame containing the serialized GenomeKit objects.
        target_cols: A dictionary mapping column names to their corresponding GenomeKit types.
    """
    genomes_exprs = [pl.col(c).struct.field("genome_str") for c in target_cols.keys()]
    genomes = (
        lf.select(
            pl.concat_list(genomes_exprs)
            .explode()
            .drop_nulls()
            .unique()
            .alias("genome_str")
        )
        .collect()["genome_str"]
        .to_list()
    )

    for genome_str in genomes:
        gk.Genome(genome_str).genes


def _validate_gkdf_metadata(metadata: dict[str, str]) -> None:
    """Validate the parquet metadata for a gk."""

    try:
        version = GkDfVersion(metadata.get("gkdf_version"))
        assert version in GkDfVersion, (
            f"Unrecognized gkdf version in Parquet metadata, expected one of {[v.value for v in GkDfVersion]}"
        )
    except ValueError:
        raise ValueError(
            "Invalid or missing gkdf_version in Parquet metadata, unable to deserialize GenomeKit objects. "
        )


def from_parquet(
    path: str, columns: list[str] | None = None, lazy: bool = False
) -> pl.DataFrame | pl.LazyFrame:
    """Deserialize a Parquet file containing GenomeKit objects into a Polars DataFrame or LazyFrame.

    Args:
        path: The file path to read the Parquet file from.
        columns: Optional list of columns to deserialize. If None, all detected
            GenomeKit columns will be deserialized.
        lazy: If True, return a LazyFrame. Otherwise, return a DataFrame.

    Returns:
        A Polars DataFrame or LazyFrame with deserialized GenomeKit objects.
    """

    metadata = pl.read_parquet_metadata(path)
    _validate_gkdf_metadata(metadata)
    target_cols = json.loads(metadata.get("target_cols"))

    lf = pl.scan_parquet(path)

    # collect unique genome strings in the file to initialize, prevents race conditions
    # on opening dganno files
    _init_gk_annotations(lf, target_cols)

    lf = lf.with_columns(
        pl.col(col)
        .map_batches(
            _map_batches_safe(REGISTRY[CURRENT_VERSION][target_cols[col]].deserializer),
            return_dtype=pl.Object,
        )
        .alias(col)
        for col in target_cols
    )

    return lf if lazy else lf.collect()
