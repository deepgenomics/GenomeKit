.. _api:

-----------------
API Documentation
-----------------

Interval
========

.. autoclass:: genome_kit.Interval
    :special-members:
    :members:

Variant
=======

.. autoclass:: genome_kit.Variant
    :special-members:
    :members:

Genome
======

.. autoclass:: genome_kit.Genome
    :special-members:
    :members:

VariantGenome
=============

.. autoclass:: genome_kit.VariantGenome
    :special-members:
    :members:

Tracks
======

GenomeDNA
---------

.. autoclass:: genome_kit.GenomeDNA
    :special-members:
    :members:

GenomeTrack
-----------

.. autoclass:: genome_kit.GenomeTrack
    :special-members:
    :members:

GenomeTrackBuilder
------------------

.. autoclass:: genome_kit.GenomeTrackBuilder
    :special-members:
    :members:

Annotations
===========

GenomeAnnotation
----------------

.. autoclass:: genome_kit.GenomeAnnotation
    :special-members:
    :members:

Gene
----

.. autoclass:: genome_kit.Gene
    :special-members:
    :members:

GeneTable
---------

.. autoclass:: genome_kit.GeneTable
    :special-members:
    :members:

Transcript
----------

.. autoclass:: genome_kit.Transcript
    :special-members:
    :members:

TranscriptTable
---------------

.. autoclass:: genome_kit.TranscriptTable
    :special-members:
    :members:

Exon
----

.. autoclass:: genome_kit.Exon
    :special-members:
    :members:

ExonTable
---------

.. autoclass:: genome_kit.ExonTable
    :special-members:
    :members:


Intron
------

.. autoclass:: genome_kit.Intron
    :special-members:
    :members:

IntronTable
-----------

.. autoclass:: genome_kit.IntronTable
    :special-members:
    :members:

Cds
---

.. autoclass:: genome_kit.Cds
    :special-members:
    :members:

CdsTable
--------

.. autoclass:: genome_kit.CdsTable
    :special-members:
    :members:

Read Alignments
===============

ReadAlignments
--------------

.. autoclass:: genome_kit.ReadAlignments
    :special-members:
    :members:

Junction
--------

.. autoclass:: genome_kit.Junction
    :special-members:
    :members:

JunctionTable
-------------

.. autoclass:: genome_kit.JunctionTable
    :special-members:
    :members:

Alignment
---------

.. autoclass:: genome_kit.Alignment
    :special-members:
    :members:

AlignmentTable
--------------

.. autoclass:: genome_kit.AlignmentTable
    :special-members:
    :members:

AlignmentMatch
--------------

.. autoclass:: genome_kit.AlignmentMatch
    :special-members:
    :members:

AlignmentMatchTable
-------------------

.. autoclass:: genome_kit.AlignmentMatchTable
    :special-members:
    :members:


Junction Read Alignments
========================

JReadAlignments
---------------

.. autoclass:: genome_kit.JReadAlignments
    :special-members:
    :members:

JunctionReadAlignmentsTable
---------------------------

.. autoclass:: genome_kit.JunctionReadAlignmentsTable
    :special-members:
    :members:

JunctionReadAlignments
----------------------

.. autoclass:: genome_kit.JunctionReadAlignments
    :special-members:
    :members:

JunctionReadAlignment
---------------------

.. autoclass:: genome_kit.JunctionReadAlignment
    :special-members:
    :members:

ReadDistributions
-----------------

.. autoclass:: genome_kit.ReadDistributions
    :special-members:
    :members:

JunctionReadDistributionTable
-----------------------------

.. autoclass:: genome_kit.JunctionReadDistributionTable
    :special-members:
    :members:

JunctionReadDistribution
------------------------

.. autoclass:: genome_kit.JunctionReadDistribution
    :special-members:
    :members:

JunctionReadCount
-----------------

.. autoclass:: genome_kit.JunctionReadCount
    :special-members:
    :members:

VCF Tables
==========

VCFTable
--------

.. autoclass:: genome_kit.VCFTable
    :special-members:
    :members:

VCFVariant
----------

.. autoclass:: genome_kit.VCFVariant
    :special-members:
    :members:

VariantTable
---------------

.. autoclass:: genome_kit.VariantTable
    :special-members:
    :members:

Data Access
===========

genome_kit.DataManager
----------------------

.. autoclass:: genome_kit.DataManager
    :special-members:
    :members:

genome_kit.DefaultDataManager
-----------------------------

.. autoclass:: genome_kit.DefaultDataManager
    :special-members:
    :members:

genome_kit.GCSDataManager
-------------------------

.. autoclass:: genome_kit.GCSDataManager
    :special-members:
    :members:

Internals
=========

genome_kit.gk_data
-------------------

.. automodule:: genome_kit.gk_data
.. autofunction:: genome_kit.gk_data.get_file
.. autofunction:: genome_kit.gk_data.upload_file
.. autofunction:: genome_kit.gk_data.wget

genome_kit.data_manager
-----------------------

.. autoclass:: genome_kit.data_manager.ProgressPercentage
    :special-members:
    :members:

genome_kit._twobitutils
-----------------------
.. autofunction:: genome_kit._twobitutils.write2bit

genome_kit._apply_variants
--------------------------
.. autofunction:: genome_kit._apply_variants.apply_variants

genome_kit._cxx_util
--------------------
.. autofunction:: genome_kit._cxx_util.mock
.. autofunction:: genome_kit._cxx_util.mock_result
.. autofunction:: genome_kit._cxx_util.mock_unreachable
.. autofunction:: genome_kit._cxx_util.strip_mock_bases
