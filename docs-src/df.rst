.. _df:

DataFrame Utilities
===================

The :py:mod:`genome_kit.df` subpackage contains utilities for working with Polars DataFrames that contain GenomeKit objects. This includes utilities for serializing DataFrames with GenomeKit objects to Parquet and deserializing them back to GenomeKit objects. This is useful when sharing tabular data sets, or when saving intermediate DataFrames to disk during data processing.

.. important::

    ``genome_kit.df`` depends on optional ``polars`` dependencies, which are not installed by default. These can be installed with the ``[df]`` extra:
    
    .. code-block:: bash

        pip install "genomekit[df]"

    The ``[df]`` extra is not included in the default ``genomekit`` installation.

    If you are running an x86 version of Python on an Apple Silicon Mac (e.g. M1 chip), this will also install the  ``polars-runtime-compat`` package, which is required to run Polars on Apple Silicon due to AVX features compatibility issues.


Quickstart
-----------
The serialization and deserialization entry points are :py:func:`~genome_kit.df.read_parquet` and :py:func:`~genome_kit.df.write_parquet`:

.. code-block:: python

    import polars as pl
    import genome_kit as gk

    genome = gk.Genome("ncbi_refseq.v110")
    df = pl.DataFrame(
        {
            "gene": [genome.genes[0], genome.genes[1]],
            "score": [0.1, 0.8],
        }
    )

    gk.write_parquet(df, "genes.parquet")
    ...
    ...
    restored_df = gk.read_parquet("genes.parquet")


.. note::
    
    The written parquet files can be read by any software that supports the parquet format, but the GenomeKit objects will only be restored when read with :py:func:`genome_kit.df.read_parquet`.
    

Supported GenomeKit Objects
---------------------------
The currently supported GenomeKit objects for serialization are:

- :py:class:`genome_kit.Genome`
- :py:class:`genome_kit.Interval`
- :py:class:`genome_kit.Transcript`
- :py:class:`genome_kit.Gene`
- :py:class:`genome_kit.Exon`
- :py:class:`genome_kit.Intron`
- :py:class:`genome_kit.CDS`
- :py:class:`genome_kit.UTR`
- :py:class:`genome_kit.Variant`

Public API
----------------
.. automodule:: genome_kit.df
    :members:
    :undoc-members:
    :show-inheritance: