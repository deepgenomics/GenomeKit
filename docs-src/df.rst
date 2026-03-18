.. _df:

GenomeKit DataFrame Utilities
=============================

The :py:mod:`genome_kit.df` subpackage contains utilities for working with Polars DataFrames that contain GenomeKit objects. This includes utilities for serializing DataFrames with GenomeKit objects to Parquet and deserializing them back to GenomeKit objects. This is useful when sharing tabular data sets, or when saving intermediate DataFrames to disk during data processing.

.. important::

    ``genome_kit.df`` depends on optional ``polars`` dependencies, which are not installed by default. This can be installed with the ``[df]`` extra:
    
    .. code-block:: bash

        mamba install "genomekit[df]"

    The ``[df]`` extra is not included in the default installation.


Quickstart
-----------
The serialization and deserialization entry points are :py:func:`~genome_kit.df.to_parquet` and :py:func:`~genome_kit.df.from_parquet`:

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

    gk.to_parquet(df, "genes.parquet")
    ...
    ...
    restored_df = gk.from_parquet("genes.parquet")


.. note::
    
    The written parquet files can be read by any software that supports the parquet format, but the GenomeKit objects will only be restored when read with :py:func:`~genome_kit.df.from_parquet`.
    

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

Public API
----------------
.. automodule:: genome_kit.df
    :members:
    :undoc-members:
    :show-inheritance: