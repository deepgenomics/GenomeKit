.. _develop:

==========
Developers
==========

Developing for GenomeKit currently requires a C++ compiler that support C++20 features. Both can be installed via conda (``gxx_linux-64`` and ```clangxx_osx-64``). For OSX, you will also need the SDK via ``xcode-select --install`` or from https://developer.apple.com/download/all/.

Setting up
----------

Clone the source tree::

    git clone git@github.com:deepgenomics/GenomeKit.git

From the ``GenomeKit`` directory, install the provided conda environment which
contains all dependencies::

    conda env create -f genomekit_dev.yml
    conda activate genomekit_dev

On Windows, you'll need to comment out the mac/linux only test dependencies from genomekit_dev.yml.

On M1 macs, you might need to set up the environment differently::

    conda create -n cxx cxx-compiler zlib fmt
    conda activate cxx
    conda install -c conda-forge -c bioconda --file a-file-with-the-deps-from-genomekit_dev-yml.txt

Build the package in development mode::

    pip install -e .

This builds the C++ extension and copies it into
your source tree (``genome_kit/_cxx.so``).
It also ensures that ``import genome_kit`` works from any directory
by linking your source tree from python's ``site-packages``.

.. note:: Windows Prerequisites

    You will need VS 2019 or newer installed. To get a compatible shell, either locate 
    and run ``vcvars64.bat``, or start the x64 Native Tools Command Prompt from the
    Start menu.
    To open VS with a preconfigured project, directly run in that command prompt::

        .vcproj\genome_kit.sln

Finally, run the all the tests::

    python -m unittest discover

You can also run examples from the ``demos`` directory.


Jetbrains CLion setup
---------------------

In the CMake settings, set the following environment variables::

    IN_CLION=1;CONDA_PREFIX=$HOME/conda/envs/genomekit_dev


Making changes
--------------

If the C/C++ code changed, you must re-run the ``develop`` command::

    pip install -e .

This includes switching branches, merging changes, or editing the C/C++ code
yourself. *Forgetting this step may lead to unpredictable behaviour.*

.. tip:: To speed up compilation on Ubuntu or Mac, install ``ccache``.

Before checking in any changes, run all tests locally::

    python -m unittest discover


Adding tests
------------

Tests are located in the ``tests`` directory, and any data they need
is located in the ``tests/data`` directory.

While developing a test, you may want to run it repeatedly, without
all other tests.
For example, to run just the ``TestInterval.test_serialize`` method in
``tests/test_interval.py`` use::

    python -m unittest tests.test_interval.TestInterval.test_serialize

C++ tests
^^^^^^^^^

To test C++ code directly, you can compile and run src/main.cpp::

    cmake -DCMAKE_BUILD_TYPE=Debug -B unittestbuild
    cmake --build unittestbuild --parallel --verbose --target main test

Debugging tests
^^^^^^^^^^^^^^^

Define envvar `GK_DEBUGBREAK` to break upon GK_CHECK failures when running
under a debugger.


Building data files
-------------------

GenomeKit relies on many pre-built files.
For example, the binary annotation ``gencode.v19.annotation.dganno``
is built from ``gencode.v19.annotation.gff3.gz``.
Reasons to re-build these files include:

* Changes to the binary file format.
* Updates to the source data.
* Changes to the processing of source data.

GenomeKit has two sets of data files:

* *Full data files* are for normal use.
  They are stored remotely in the GenomeKit store
  and pulled to the user's local file system on-demand.

* *Test data files* are for testing.
  They are tiny excerpts of the full files, small enough
  to check in to source control, fast enough to run in
  continuous integration testing.
  They are stored in the source tree under ``tests/data``.

The ``genome_kit`` module's ``build`` command can be used to build full
Appris/MANE data files, and Appris/MANE/dganno/2bit test data files.

For a full set of options, run::

    python -m genome_kit build --help


Building full data files
^^^^^^^^^^^^^^^^^^^^^^^^

For instructions on how to build annotation (dganno) files and assembly
(2bit) files, see `Genomes <genomes.html>`_.

Full-sized data files reside in a local user directory reserved
for GenomeKit, downloaded from the data store on-demand.

.. note:: See the API Documentation for instructions on how to build
    `data tracks <api.html#genometrackbuilder>`_,
    `read alignments <api.html#genome_kit.ReadAlignments.build_ralign>`_,
    `read distributions <api.html#genome_kit.ReadDistributions.build_rdist>`_,
    `junction read alignments <api.html#genome_kit.JReadAlignments.build_jralign>`_,
    and `VCF tables <api.html#genome_kit.VCFTable.build_vcfbin>`_.


Building test data files
^^^^^^^^^^^^^^^^^^^^^^^^

Test data files reside in the source tree under ``tests/data``.
To build them, you must have registered your source tree in
develop mode::

    pip install -e .

Now that your source tree is the default `genome_kit` import,
the ``build`` subcommand will be able to find
your test data directory.

To build test annotation, 2bit, Appris, and MANE files, use `--test-<type>`
flags on the ``build`` subcommand::

    python -m genome_kit build --test-anno --test-2bit --test-appris --test-mane


Releasing GenomeKit
-------------------

The `GenomeKit repo <https://github.com/deepgenomics/GenomeKit>`__ uses
the `Release Please bot <https://github.com/googleapis/release-please>`__
to create Github releases based on PRs. When the bot creates a PR, you can
merge it to create a release.

Once a Github release is created, a PR will automatically be created in
the `GenomeKit conda-forge feedstock repo <https://github.com/conda-forge/genomekit-feedstock>`__
by regro-cf-autotick-bot. Once that PR is merged, conda-forge's CI
pipeline is kicked off and the new version of GenomeKit is built and published
to conda-forge.
