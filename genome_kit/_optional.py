from __future__ import annotations


def require_polars():
    """Import Polars if available, otherwise fail gracefully."""
    try:
        import polars as pl
    except ModuleNotFoundError as e:
        raise ImportError(
            "Optional dependency 'polars' is required for this functionality. Please install with `pip install genomekit[df]`."
        ) from e

    return pl
