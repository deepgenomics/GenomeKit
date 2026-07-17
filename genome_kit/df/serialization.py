from __future__ import annotations

import functools
import itertools
import json
import warnings
from collections.abc import Callable
from inspect import signature
from pathlib import Path
from typing import TYPE_CHECKING, Any, overload, TypeVar, TypeAlias

if TYPE_CHECKING:
    # import libraries for static type checkers
    import pandas as pd
    import polars as pl

    # supported dataframe types for serialization/deserialization
    SupportedTabular: TypeAlias = pl.DataFrame | pl.LazyFrame | pd.DataFrame
    DF = TypeVar("DF", pl.DataFrame, pl.LazyFrame, pd.DataFrame)

import genome_kit as gk
from genome_kit._optional import require_pandas, require_polars

from .gk_structs import CURRENT_VERSION, CellType, ColumnInfo, GkDfType, GkDfVersion, identify_struct
from .registry import GK_TO_GKDF_TYPE, get_registry


def _map_batches_safe(fn: Callable) -> Callable:
    """Wrap a user defined function (UDF) and run safely with polars map_batches.

    Polars has a bug in map_batches that incorrectly forwards the return_dtype argument
    to the UDF. See https://github.com/pola-rs/polars/issues/24840.
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
    """Infer columns containing GenomeKit objects and their shape (list or scalar).

    Uses the first `infer_schema_length` rows for inference.
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
        head_types = {type(v) for v in vals}

        if isinstance(first, list):
            if head_types != {list}:
                raise ValueError(
                    f"Column {col} contains mixed data types: {list(itertools.islice(head_types, 3))}.\n"
                    "Please ensure all cells are the same type before serialization."
                )
            cell_type = CellType.LIST
            col_types = {type(item) for v in vals for item in v if item is not None}
        else:
            cell_type = CellType.SCALAR
            col_types = set(vals.map_elements(type, return_dtype=pl.Object))

        if len(col_types) != 1:
            raise ValueError(
                f"Column {col} contains mixed data types: {list(itertools.islice(col_types, 3))}.\n"
                "Please ensure all cells are the same type before serialization."
            )

        col_type = GK_TO_GKDF_TYPE.get(col_types.pop(), None)

        if col_type is None:
            # column is not a genomekit type, so no serialization needed
            continue

        target_cols[col] = ColumnInfo(cell_type=cell_type, gkdf_type=col_type)

    return target_cols


def _list_serializer(
    serializer: Callable[[pl.Series], pl.Series], return_dtype: Any
) -> Callable[[pl.Series], pl.Series]:
    """Convert a serializer to accept a series of lists of objects.
    
    Default serializers accept a series of single objects.
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
    """
    pl = require_polars()

    def genome_str_field(col_info: dict) -> str:
        gkdf_type = col_info["gkdf_type"]
        if gkdf_type == GkDfType.GENOME:
            return "genome_name"
        elif gkdf_type in (GkDfType.INTERVAL, GkDfType.VARIANT):
            return "refg"
        else:
            return "anno"

    anno_strong_refs = []

    # extract genome_str field from every column
    genomes_exprs = []
    genomes_list_exprs = []

    for c in target_cols.keys():
        genome_field = genome_str_field(target_cols[c])
        if target_cols[c]["cell_type"] == CellType.SCALAR:
            genomes_exprs.append(pl.col(c).struct.field(genome_field))
        else:
            genomes_list_exprs.append(pl.col(c).explode().struct.field(genome_field))

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
                pl.concat(genomes_list_exprs)
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
        try:
            anno_strong_refs.append(genome.genes)
        except ValueError:
            # reference genomes don't have annotations
            continue

    return anno_strong_refs


def _validate_gkdf_metadata(metadata: dict[str, str]) -> None:
    # gkdf version
    metadata_version = metadata.get("gkdf_version")
    version = GkDfVersion(metadata_version) if metadata_version is not None else None
    if version != CURRENT_VERSION:
        raise ValueError(
            f"Invalid or missing gkdf_version in Parquet metadata, unable to deserialize GenomeKit objects. "
            f"Expected GkDfVersion {CURRENT_VERSION}, but found {version}."
        )

    # target cols
    if metadata.get("target_cols") is None:
        raise ValueError(
            "Missing target_cols in Parquet metadata, unable to deserialize GenomeKit objects."
        )

    # gk version
    gk_version = metadata.get("gk_version")
    if gk_version is None:
        raise ValueError("Missing gk_version in Parquet metadata.")
    elif gk_version != gk.__version__:
        warnings.warn(
            f"Parquet file was written with GenomeKit version {gk_version}, but current version is {gk.__version__}. "
            "Deserializing GenomeKit objects may not be consistent across versions."
        )


def _list_deserializer(
    deserializer: Callable[[pl.Series], pl.Series],
) -> Callable[[pl.Series], pl.Series]:
    """Convert a deserializer to accept a series of lists of objects.

    Default deserializers accept a series of single objects.
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
    """Deserialize specified columns containing GenomeKit objects.

    target_cols is a dictionary representation of the ColumnInfo dataclass.
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


def _convert_pandas_to_polars(df: pd.DataFrame) -> pl.LazyFrame:
    """Convert a pandas DataFrame to a Polars LazyFrame.

    Uses an intermediate representation to remove dependency on pyarrow for conversion.
    """
    pl = require_polars()

    # pandas allows duplicate column names
    if any(len(df[col].shape) > 1 for col in df.columns):
        raise ValueError(
            "Input DataFrame contains duplicated column names. "
            "Unique column names are required for serialization."
        )

    lf = pl.LazyFrame(df.to_dict(orient="list"), strict=False)
    # fill np.nan with nulls for consistent "None" values in polars
    # ONLY applies to float columns, np.nan in object (GenomeKit) columns will remain
    lf = lf.fill_nan(None)

    return lf


def _convert_to_polars_lf(df: SupportedTabular) -> pl.LazyFrame:
    pl = require_polars()

    if isinstance(df, pl.DataFrame):
        return df.lazy()
    elif isinstance(df, pl.LazyFrame):
        return df
    # passed object is not a polars DataFrame/LazyFrame, check module and import
    if type(df).__module__.startswith("pandas"):
        pd = require_pandas()
        if isinstance(df, pd.DataFrame):
            return _convert_pandas_to_polars(df)

    raise TypeError(
        f"Unsupported DataFrame type {type(df)}. Please provide a Polars DataFrame or LazyFrame, or a pandas DataFrame."
    )


def write_parquet(
    df: SupportedTabular,
    path: str | Path,
    infer_schema_length: int = 100,
) -> None:
    """Serialize a DataFrame or LazyFrame with GenomeKit objects to a Parquet file.

    Args:
        df: A Polars DataFrame or LazyFrame or pandas DataFrame with columns containing GenomeKit objects.
        path: The file path to write the Parquet file to.
        infer_schema_length: The number of rows to use for schema inference when writing the Parquet file.
    """
    pl = require_polars()

    path = Path(path)
    # convert input to a polars LazyFrame for processing.
    df = _convert_to_polars_lf(df)

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


def _process_genomekit_parquet(path: Path, deserialize_gk_objects: bool = True) -> pl.LazyFrame:
    pl = require_polars()
    metadata = pl.read_parquet_metadata(path)
    _validate_gkdf_metadata(metadata)
    target_cols = json.loads(metadata.get("target_cols"))

    lf = pl.scan_parquet(path)

    if deserialize_gk_objects:
        # collect unique genome strings in the file and initialize, prevents race conditions
        # on opening dganno files in concurrent scenarios.
        # genomes returned in dummy variable to keep weak reference alive for deserialization
        _ = _init_gk_annotations(lf, target_cols)

        lf = _deserialize_gk_cols(lf, target_cols)

    return lf


def _convert_to_output_format(lf: pl.LazyFrame, astype: type[DF]) -> DF:
    pl = require_polars()

    if astype is pl.DataFrame:
        return lf.collect()
    elif astype is pl.LazyFrame:
        return lf
    elif astype.__module__.startswith("pandas"):
        pd = require_pandas()
        if astype is pd.DataFrame:
            return lf.collect().to_pandas(use_pyarrow_extension_array=False)

    raise TypeError(
        f"Unsupported astype {astype}. Please provide pl.DataFrame, pl.LazyFrame, or pd.DataFrame."
    )


@overload
def read_parquet(path: str | Path) -> pl.DataFrame: ...

@overload
def read_parquet(path: str | Path, astype: type[DF]) -> DF: ...


def read_parquet(path: str | Path, astype: type[DF] | None = None, deserialize_gk_objects: bool = True) -> DF:
    """Deserialize a Parquet file containing GenomeKit objects into a tabular data format.

    The type of the returned object is determined by the `astype` argument.

    Args:
        path: The file path to read the Parquet file from.
        astype: The data type of tabular data to return. Defaults to a Polars DataFrame.

    Returns:
        A tabular data format with the deserialized GenomeKit objects.
    """
    pl = require_polars()

    path = Path(path)
    lf = _process_genomekit_parquet(path, deserialize_gk_objects)

    return _convert_to_output_format(lf, astype or pl.DataFrame)



