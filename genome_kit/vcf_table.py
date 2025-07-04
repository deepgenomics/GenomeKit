# Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
from __future__ import annotations

import contextlib
import os
import subprocess
import tempfile
import shutil
from typing import Optional, Union
from typing_extensions import Literal

import numpy as np

from . import _cxx
from ._cxx_util import mock, mock_result, mock_unreachable, strip_mock_bases
from .genome import Genome
from .variant import Variant


def as_vcf_chrom(chrom):
    return chrom.replace("chrM", "chrMT").replace("chr", "")


@contextlib.contextmanager
def _dump_fasta(genome, chrom_mapper):
    """Dumps a FASTA containing the DNA sequence of refg, indexes it, and returns a context
    that can be treated as the path to the fasta. The temporaries are deleted once the
    context is closed."""
    samtools = shutil.which("samtools")
    if not samtools:
        raise RuntimeError("samtools must be installed in order to normalize VCFs")  # pragma: no cover

    with tempfile.NamedTemporaryFile(mode="w", prefix="%s_" % genome.refg, suffix=".fasta") as f:
        for chrom in genome.chromosomes:
            try:
                chrom_size = genome.chromosome_size(
                    chrom)  # Currently may error out or may miss chromosomes for non-human species
            except ValueError:
                pass  # pragma: no cover
            if not chrom_size:
                continue

            chunk_size = 10000
            f.write(">%s\n" % chrom_mapper(chrom))
            for i in range(0, chrom_size, chunk_size):
                interval = genome.interval(chrom, "+", i, min(i + chunk_size, chrom_size))
                f.write(genome.dna(interval))
                f.write("\n")

        f.flush()

        try:
            subprocess.check_output([samtools, "faidx", f.name], stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            raise RuntimeError("{}: {}".format(e, e.output))  # pragma: no cover

        try:
            yield f.name
        finally:
            os.unlink(f.name + ".fai")  # Also clean up index


@contextlib.contextmanager
def _dump_vcfbuild_script(vcfbinpath, refg, *args, **kwargs):
    """Dumps a python script that will build a specific vcfbin file from stdin.
    Used to redirect output of `bcftools norm` directly to VCFTable.build_vcfbin.
    (It is not easy to redirect to stdin of the current process, yet redirecting
    is important. Creating a temporary VCF file is too large or too slow (gzip)."""
    with tempfile.NamedTemporaryFile(mode="w+", prefix="genomekit_from_vcf_", suffix=".py") as script:
        script.write(f"""
import sys
from genome_kit import Genome, VCFTable, Interval
vcfbinpath = '{vcfbinpath}'
refg = Genome('{refg.config}')
VCFTable.build_vcfbin(vcfbinpath,
                      sys.stdin,
                      refg,
                      {",".join(repr(x) for x in args)}{"," if args else ""}
                      {",".join(f"{k}={repr(v)}" for (k, v) in kwargs.items())})""")
        script.flush()
        yield script.name


@strip_mock_bases
@_cxx.register
class VCFVariant(_cxx.VCFVariant, Variant):
    """A VCF variant. Scalar INFO/FORMAT fields can be accessed as a dynamic
    attribute, eg. ``variant.AD``.

    Bases: :py:class:`~genome_kit.Variant`
    """
    __slots__ = ()  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    # noinspection PyMissingConstructor
    @mock
    def __init__(self):  # pragma: no cover
        pass  # Stub to prevent Interval.__init__ from being accessible to docs

    def as_variant_string(self):
        """Returns a string representation of the variant.

        See Also
        --------
        :py:meth:`~genome_kit.Variant.from_string`

        Returns
        -------
        :py:class:`str`
            The variant in the format ``"chromosome:start:ref:alt"``, where
            start is in DNA1 coordinates.
        """
        return '{}:{}:{}:{}'.format(self.chromosome, self.start + 1, self.ref, self.alt)

    __str__ = as_variant_string

    def __repr__(self):
        return '<VCFVariant {}:{}:{}:{}>'.format(self.sys, self.start + 1, self.ref, self.alt)


@_cxx.register
class VCFTable(_cxx.VCFTable):
    """Open a ``.vcfbin`` file for querying.

    See the VCF section of :ref:`quickstart` for examples of how to use this class,
    or see ``demos/query_vcf.py``.

    To convert from ``.vcf[.gz]`` see :py:meth:`~genome_kit.VCFTable.build_vcfbin`.

    Parameters
    ----------
    infile : :py:class:`str`
        The path to the ``.vcfbin`` file.
    """

    __slots__ = ("_mask_indices",)  # <--- NEW SLOTS OK
    # def __new__(self):  pass   # <--- DO NOT OVERRIDE BASE CLASS IMPLEMENTATION
    # def __del__(self):  pass   # <--- DO NOT IMPLEMENT

    def __init__(self, infile):
        self._mask_indices = None
        _cxx.VCFTable.__init__(self, infile)

    @staticmethod
    def from_vcf(vcfpath: str,
                 reference_genome: Genome,
                 *args,
                 normalize: Union[bool, Literal["EBI", "UCSC"]] = False,
                 cache: bool = True,
                 vcfbinpath: Optional[str] = None,
                 **kwargs) -> VCFTable:
        """
        Returns a :py:class:`~genome_kit.VCFTable` object opened to a VCF file.

        The function internally converts the VCF to a binary format and then returns a `VCFTable` object
        for fast access to the variants and fields within. By default the binary file (`.vcf.bin`) will be
        written adjacent to the input file.

        It is intended to be used as follows::

            def open_clinvar(vcfpath):
                return VCFTable.from_vcf(vcfpath, "hg19", info_ids=["AF_EXAC", "CLNSIG"])

            vcf = open_clinvar("clinvar_20180128.vcf.gz")
            variants = vcf.find_within(brca1_interval)
            ...

        If the input VCF contains multi-allelic sites, then use `normalize` to automatically
        split them on conversion. This step requires `samtools` and `bcftools` to be installed.

        See Also
        --------
        :py:meth:`~genome_kit.VCFTable.build_vcfbin` for all the possible parameters.

        Parameters
        ----------
        vcfpath
            Path to input `.vcf or `.vcf.gz` file.

        reference_genome
            The reference genome (e.g. `"hg19"`). Forwarded to :py:meth:`~genome_kit.VCFTable.build_vcfbin`.

        args
            Parameters forwarded to :py:meth:`~genome_kit.VCFTable.build_vcfbin`.

        normalize
            Whether to normalize the VCF during conversion (if applicable).
            If `True`, `bcftools` is used for normalization of all the EBI
            (Ensembl) chromosomes; alternatively, `EBI` or `UCSC` can be used
            to specify the naming scheme used by the VCF (eg. GTEx uses the
            UCSC format).
            Unsupported on Windows due to bcftools and samtools lack of support.

        cache
            Whether to re-use the binary or to re-generate it from `vcfpath`.

        vcfbinpath
            Path to the intermediate binary. Default is adjacent to `vcfpath`.

        kwargs
            Parameters forwarded to :py:meth:`~genome_kit.VCFTable.build_vcfbin`.

        Returns
        -------
        :py:class:`~genome_kit.VCFTable`
            The opened VCF file.
        """
        # If user doesn't specify path for binary output file, write it adjacent to the vcf
        if vcfbinpath is None:
            vcfbinpath = vcfpath.replace(".gz", "") + ".bin"  # pragma: no cover

        if cache and os.path.exists(vcfbinpath):
            if os.path.getmtime(vcfpath) > os.path.getmtime(vcfbinpath):
                print("Warning: %s older than %s in VCFTable.from_vcf." % (vcfbinpath, vcfpath))  # pragma: no cover
        else:
            if normalize:
                # Check requirements
                bcftools = shutil.which("bcftools")
                if not bcftools:
                    raise RuntimeError("bcftools must be installed to normalize VCFs")  # pragma: no cover

                if normalize not in [True, "EBI", "UCSC"]:
                    raise ValueError('normalize must be in [True, "EBI", "UCSC"]')
                chrom_mapper = as_vcf_chrom if normalize != "UCSC" else lambda x: x
                regions = ",".join(map(chrom_mapper, reference_genome.chromosomes))

                # Start one subprocess to normalize the VCF using the genome's fasta, and pipe the
                # stdout of that into another subprocess that directly builds the vcfbin.
                with _dump_fasta(reference_genome, chrom_mapper) as fasta:
                    with _dump_vcfbuild_script(vcfbinpath, reference_genome, *args, **kwargs) as script:
                        norm = subprocess.Popen(
                            [bcftools, "norm", "-r", regions, "-m", "-", "-c", "wx", "-f", fasta, vcfpath],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
                        vcfbin = subprocess.Popen(["python", script],
                                                  stdin=norm.stdout,
                                                  stdout=None,
                                                  stderr=subprocess.PIPE)
                        norm.stdout.close()
                        _, err = vcfbin.communicate()
                        if vcfbin.returncode != 0:
                            raise RuntimeError(err)
                        if norm.wait() != 0:
                            raise RuntimeError("bcftools failed: {}".format(norm.stderr.read()))
            else:
                VCFTable.build_vcfbin(
                    vcfbinpath,
                    vcfpath,
                    reference_genome,
                    *args,
                    **kwargs)

        # Open the binary file
        return VCFTable(vcfbinpath)

    def find_5p_aligned(self, interval):
        """Returns all variants that have 5' end aligned with 5' end of `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant.end5 == interval.end5]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that have 5' end aligned to the 5' end of `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_5p_aligned(self, interval))

    def find_3p_aligned(self, interval):
        """Returns all variants that have 3' end aligned with 3' end of `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant.end3 == interval.end3]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that have 3' end aligned to the 3' end of `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_3p_aligned(self, interval))

    def find_5p_within(self, interval):
        """Returns all variants that have 5' end within `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant.end5 in interval]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that have 5' end falling entirely within `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_5p_within(self, interval))

    def find_3p_within(self, interval):
        """Returns all variants that have 3' end within `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant.end3 in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that have 3' end falling entirely within `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_3p_within(self, interval))

    def find_within(self, interval):
        """Returns all variants that fall entirely within `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant in interval]

        The results will also be sorted by 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that fall entirely within `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_within(self, interval))

    def find_overlapping(self, interval):
        """Returns all variants that overlap `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant.overlaps(interval)]

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that overlap `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_overlapping(self, interval))

    def find_exact(self, interval):
        """Returns all variants that span exactly `interval`.

        This function is a faster version of the brute-force filter::

            [variant for variant in vcf if variant.interval == interval]

        The results will also be sorted by 5' and then 3'.

        Parameters
        ----------
        interval : :py:class:`~genome_kit.Interval`
            The query interval.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
            All variants that span exactly `interval`.
        """
        return self._masked_results(_cxx.VCFTable.find_exact(self, interval))

    @mock
    @property
    def stranded(self):  # pragma: no cover
        """If `True` then the strand is significant when calling the `find_x` methods.

        Returns
        -------
        :py:class:`bool`
            Whether this table can contain negative stranded intervals.
        """
        return mock_result(bool)

    def __getitem__(self, index):
        """Access to all variants.

        Allows iteration over all variants::

            for variant in variants:
                ...

        The results will be ordered as in the source file.

        Parameters
        ----------
        index : :py:class:`int`
            The integer index of the requested variant.

        Returns
        -------
        :py:class:`~genome_kit.VCFVariant`
           The variant object identified by the given index.

        Notes
        -----
        The VCFVariant will be invalidated after this table has been closed via :py:meth:`~.close` or
        a ``with`` statement.
        """
        if self._mask_indices is not None:
            index = self._mask_indices[index]
        return _cxx.VCFTable.__getitem__(self, index)

    def __len__(self):
        return _cxx.VCFTable.__len__(self) if self._mask_indices is None else len(self._mask_indices)

    def where(self, mask):
        """Filter variants by numpy mask.

        Fast extraction of table elements based on a mask::

            # Intended to be used like this.
            variants = vcf.where(mask)

            # It is faster but equivalent to this.
            variants = [vcf[i] for i in np.where(mask)[0]]

        Parameters
        ----------
        mask : :py:class:`ndarray`
            A boolean mask the same size as the table.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
           A list of variants selected by the mask.
        """
        if self._mask_indices is not None:
            indices = self._mask_indices[mask]
            mask = np.zeros(_cxx.VCFTable.__len__(self), np.bool_)
            mask[indices] = True
        return _cxx.VCFTable.where(self, mask)

    def masked(self, mask):
        """Return a new VCFTable which is filtered by a numpy mask

        For example, to do a random downsample::

            downsampled = vcf.masked(np.random.binomial(1, 0.1, len(vcf)).astype(bool))

        Or to subset by only ABCA4 genes::

            abca4_vcf = vcf.masked(vcf.index_of(v) for v in vcf.find_overlapping(gene))

        Notes
        -----
        The masked VCFTable does not share resources with the original, so :meth:`~genome_kit.VCFTable.close`
        should be called on both instances when required.

        Parameters
        ----------
        mask : :class:`ndarray` of :class:`bool` or :class:`ndarray` of :class:`int`
            A boolean mask the same size as the table or a list of indices.

        Returns
        -------
        :py:class:`list` of :py:class:`~genome_kit.VCFVariant`
           A list of variants selected by the mask.
        """
        result = VCFTable(self.filename)
        if np.dtype(type(mask[0])) is np.dtype(bool):
            result._mask_indices = np.nonzero(mask)[0] if self._mask_indices is None else self._mask_indices[mask]
        else:
            result._mask_indices = np.array(mask) if self._mask_indices is None else np.intersect1d(
                mask, self._mask_indices, assume_unique=True)
        return result

    def _masked_results(self, results):
        if self._mask_indices is None:
            return results
        indices = [_cxx.VCFTable.index_of(self, x) for x in results]
        indices = np.intersect1d(indices, self._mask_indices, assume_unique=True)
        mask = np.zeros(_cxx.VCFTable.__len__(self), np.bool_)
        mask[indices] = True
        return _cxx.VCFTable.where(self, mask)

    # values for GT (as per cyvcf2)
    GT_HOMOZYGOUS_REF = np.int8(0)
    """Constant indicating a ``0/0`` or ``0|0`` in GT format value."""
    GT_HETEROZYGOUS_UNPHASED = np.int8(1)
    """Constant indicating a ``0/1`` in GT format value."""
    GT_HOMOZYGOUS_ALT = np.int8(2)
    """Constant indicating a ``1/1`` or ``1|1`` in GT format value."""
    GT_UNKNOWN = np.int8(3)
    """Constant indicating a ``?/?`` or ``?|?`` in GT format value."""
    GT_HETEROZYGOUS_PHASED_0_1 = np.int8(4)
    """Constant indicating a ``0|1`` in GT format value."""
    GT_HETEROZYGOUS_PHASED_1_0 = np.int8(5)
    """Constant indicating a ``1|0`` in GT format value."""

    # values for SVTYPE
    SVTYPE_NA = 0
    """Constant indicating the SVTYPE was missing."""
    SVTYPE_DEL = 1
    SVTYPE_INS = 2
    SVTYPE_DUP = 3
    SVTYPE_INV = 4
    SVTYPE_CNV = 5
    SVTYPE_BND = 6

    @mock
    @property
    def num_samples(self):  # pragma: no cover
        """The number of samples represented in the VCF.

        Returns
        -------
        :py:class:`int`
            Number of samples represented in the VCF file, for which FORMAT data was provided.
        """
        return mock_result(int)

    @mock
    @property
    def sample_names(self):  # pragma: no cover
        """The names of the samples represented in the VCF if it was encoded with format IDs.

        Returns
        -------
        :class:`list` of :class:`str`
            A list of the sample IDs in the same order as the original VCF. The indices correspond to the format values.
        """
        return mock_result(list)

    @mock
    @property
    def info_ids(self):  # pragma: no cover
        """The INFO fields available for :py:meth:`~.info` that were encoded in :meth:`~.build_vcfbin` or
        :meth:`~.from_vcf`.

        Returns
        -------
        :class:`list` of :class:`str`
            A list of INFO ids that can be passed into :py:meth:`~.info`.
        """
        return mock_result(list)

    def info(self, info_id):
        """The full per-variant data array associated with an INFO field. String INFO fields
        are accessed by :class:`~.VCFVariant` attribute, eg. ``variant.CIGAR``.

        Parameters
        ----------
        info_id : :py:class:`str`
            The INFO ID from the VCF, such as `AD` or `DP`.

        Returns
        -------
        :py:class:`~numpy.ndarray`
            A complete array of INFO data values associated with `info_id`.

            If there is one value per variant, the results is 1-dimensional::

                INFO
                DP=1
                DP=4

                [1 4]

            If there are multiple values per variant, the results is 2-dimensional::

                INFO
                AD=0,1
                AD=2,2

                [[0 1]
                 [2 2]]

        Raises
        ------
        KeyError
            Raised if `info_id` is unrecognized, because it was not in the original VCF
            or included by :py:meth:`~genome_kit.VCFTable.build_vcfbin`
        TypeError
            Raised if `info_id` refers to a String valued column. String values are
            accessed by :class:`~.VCFVariant` attribute.

        """
        results = _cxx.VCFTable.info(self, info_id)
        if self._mask_indices is not None:
            results = results[self._mask_indices]
        return results

    @mock
    @property
    def format_ids(self):  # pragma: no cover
        """The INFO fields available for :py:meth:`~.format` that were encoded in :meth:`~.build_vcfbin` or
        :meth:`~.from_vcf`.

        Returns
        -------
        :class:`list` of :class:`str`
            A list of FORMAT ids that can be passed into :py:meth:`~.format`.
        """
        return mock_result(list)

    def format(self, format_id):
        """The full per-variant per-sample data matrix associated with a FORMAT field.

        Parameters
        ----------
        format_id : :py:class:`str`
            The FORMAT ID from the VCF, such as `AD` or `GT`.

        Returns
        -------
        :py:class:`~numpy.ndarray`
            A complete array of FORMAT data values associated with `format_id`.

            If there is one value per variant and per sample, the results is 2-dimensional::

                sample1 sample2 sample3
                1       2       3
                4       5       6

                [[1 2 3]
                 [4 5 6]]

            If there are multiple values per variant and per sample, the results is 3-dimensional::

                sample1 sample2 sample3
                1,2     3,4     5,6
                7,8     9,10    11,12

                [[[1 2] [3 4]  [5  6 ]]
                 [[7 8] [9 10] [11 12]]]

        Raises
        ------
        KeyError
            Raised if `format_id` is unrecognized, either because it was not in the original VCF
            or included by :py:meth:`~genome_kit.VCFTable.build_vcfbin`.

        """
        results = _cxx.VCFTable.format(self, format_id)
        if self._mask_indices is not None:
            results = results[self._mask_indices]
        return results

    @mock
    @staticmethod
    def build_vcfbin(outfile,
                     infile,
                     reference_genome,
                     info_ids=None,
                     fmt_ids=None,
                     validate=True,
                     exclude=None,
                     allow=None,
                     ancestral="error"):  # pragma: no cover
        """Builds a binary ``.vcfbin`` file from a ``.vcf[.gz]``.

        The resulting file can be opened by :py:class:`~genome_kit.VCFTable`.

        **VCF normalization.** The input VCF should be normalized, and multi-allelic sites must be split.
        Consider running `infile` through ``bcftools norm -m - -f``. Examples::

            REF    ALT
            CT     CAT     # FAIL (not normalized)
                   A       # OK
            .      A       # OK
            -      A       # OK
            C      CA      # OK
            C      A,G     # FAIL (multi-allelic)

        **Extracting INFO/FORMAT columns.** By default, INFO and FORMAT columns such as DP and AD are omitted.
        Use `info_ids` and `fmt_ids` arguments to specify which columns you want included::

            ["GT", "AD"]                        # Keep GT and AD (as int8 and int32 respectively)
            {"GT" : None, "AD" : None}          # Same
            {"GT" : None, "AD" : -1}            # Same, but store -1 for any missing AD value
            {"GT" : None, "AD" : np.int16(0) }  # Keep GT and AD (as int8 and int16 respectively)

        The GT (genotype) FORMAT and SVTYPE (structural variant type) INFO columns are special. They default to
        :py:class:`~numpy.int8`, and the string values are stored as integer constants:

        ***GT***

        * ``0/0`` and ``0|0`` are stored as :const:`VCFTable.GT_HOMOZYGOUS_REF` (0)
        * ``0/1`` is stored as :const:`VCFTable.GT_HETEROZYGOUS_UNPHASED` (1)
        * ``1/1`` and ``1|1`` are stored as :const:`VCFTable.GT_HOMOZYGOUS_ALT` (2)
        * ``?/?`` and ``?|?`` are stored as :const:`VCFTable.GT_UNKNOWN` (3)
        * ``0|1`` is stored as :const:`VCFTable.GT_HETEROZYGOUS_PHASED_0_1` (4)
        * ``1|0`` is stored as :const:`VCFTable.GT_HETEROZYGOUS_PHASED_1_0` (5)

        ***SVTYPE***

        * A missing SVTYPE is represented as :py:const:`VCFTable.SVTYPE_NA` (0)
        * ``DEL`` is stored as :const:`VCFTable.SVTYPE_DEL` (1)
        * ``INS`` is stored as :const:`VCFTable.SVTYPE_INS` (2)
        * ``DUP`` is stored as :const:`VCFTable.SVTYPE_DUP` (3)
        * ``INV`` is stored as :const:`VCFTable.SVTYPE_INV` (4)
        * ``CNV`` is stored as :const:`VCFTable.SVTYPE_CNV` (5)
        * ``BND`` is stored as :const:`VCFTable.SVTYPE_BND` (6)

        The ID INFO column is stored as an INFO ID. If it is required (e.g., you need to find a MATEID pair), pass
        `"ID"` into `info_ids`.

        Missing values are also treated differently: GT values must be specified as a diploid.
        If a single missing allele (e.g. ``0/.``, ``./0``, ``1/.``, or ``./1``) should be encoded as a *possible*
        :py:const:`VCFTable.GT_HETEROZYGOUS_UNPHASED` instead of :py:const:`VCFTable.GT_UNKNOWN`, the missing value should
        be set to :py:const:`VCFTable.GT_HETEROZYGOUS_UNPHASED`. ``./.`` will still be interpreted as
        :py:const:`VCFTable.GT_UNKNOWN`::

            {"GT" : VCFTable.GT_HETEROZYGOUS_UNPHASED}

        For phased values, the default value is ignored. For example, ``.|1`` will be encoded as
        :py:const:`VCFTable.GT_UNKNOWN`.

        All other columns default to :py:class:`~numpy.int32`, :py:class:`~numpy.float32`, or :py:class:`str` as defined
        by the VCF Type header (can override bitwidth for performance or to save disk space.)
        For integer/float columns, only constant value dimensions are currently supported.
        By default, a missing value is substituted with 0, 0.0, or "". If missing integer values must be
        distinguishable, consider using ``-(2**31)`` to ``(2**31+7)`` as these are reserved by VCF.

        **String-valued columns** String values can be accessed by attribute (``variant.CIGAR``).
        Unlike numeric columns, a string column cannot be accessed as a numpy array.

        Caveats: Multi-dimensional values are not split, and ASCII-encoded VCF values,
        eg. ``%3A`` for ``:``, will not be decoded. ``.`` are not decoded as empty strings
        unless it is the sole value.

        **Value multiplicity.** Each INFO/FORMAT field has an expected number of values per sample,
        following the VCF 4.2 Specification sections 1.4.2 and 1.6.2. For example::

            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">

        Here, a FORMAT value of ``0/1:5,8`` for a single sample would be valid because we expect one
        GT value per (``Number=1``) and two AD values (``Number=R``).

        FORMAT values with multiple samples or INFO/FORMAT values that are multidimensional can only
        be accessed via :meth:`~.info` or :meth:`~format` and cannot be accessed by :class:`~.VCFVariant`
        attributes.

        Beware VCF headers that incorrectly specify ``Number=.`` instead of ``R``, as their
        multi-allelic variants will be split incorrectly by ``bcftools``.


        Parameters
        ----------
        outfile : :py:class:`str`
            The binary output file. This file can then be opened by :py:class:`~genome_kit.VCFTable`.

        infile : :py:class:`str`
            The path to the source ``.vcf[.gz]`` file.

        reference_genome : :py:class:`~genome_kit.Genome`
            The reference genome that `infile` uses as ``#reference`` field.
            The parameter is not validated against the VCF, so it is up to the caller to ensure it is correct.
            For example, a ``#reference=b37.fa`` or ``=NCBI37`` is a match to ``"hg19"``.

        info_ids: :py:class:`list` of :py:class:`str` | :py:class:`dict` of :py:class:`str` to :py:class:`object`, \
                optional
            Controls how INFO columns are preserved in `outfile`. See above.

        fmt_ids: :py:class:`list` of :py:class:`str` | :py:class:`dict` of :py:class:`str` to :py:class:`object`, \
                optional
            Controls how FORMAT columns are preserved in `outfile`. See above.

        validate: :py:class:`bool`, optional
            Whether to validate `REF` sequences against the reference genome. Default is True.

        exclude : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Variants within these intervals will be excluded.

        allow : :py:class:`list` of :py:class:`~genome_kit.Interval`
            Optional. Only variants within these intervals will be included,
            so long as they are not excluded.

        ancestral: :py:class:`str`
            Behaviour when encountering ancestral alleles. Valid values are "error", "warn", and "exclude". Defaults
            to "error".

        Raises
        ------
        ValueError
            Under the following conditions:


                * The VCF is malformed.
                * A row in `infile` does not contain every FORMAT ID specified by `fmt_ids`.
                * An INFO/FORMAT ID is unsupported.
                * A `REF` or `ALT` sequence contained a VCF non-standard nucleotide code.
                * A `REF` sequence failed to validate.
                * An INFO/FORMAT field type is incompatible with the VCF Type header definition.

        """
        mock_unreachable()

    @mock
    @staticmethod
    def vcfbin_version():  # pragma: no cover
        """The `vcfbin` file format version that this build supports.

        Attempting to open an `vcfbin` file of a mismatched version will raise an :py:exc:`IOError`.

        Returns
        -------
        :py:class:`int`
           The version number.
        """
        return mock_result(int)

    def sequence_variations(self, interval):
        """Find variants that could affect the given `interval`'s DNA sequence.
        When specified with an `anchor`, deletion variants will increase the
        interval's span; this can be overestimated (due to overlapping
        deletions), but the spurious variants will be filtered during DNA
        sequence retrieval.

        The results will also be sorted by 5'.

        Parameters
        ----------
        interval : :class:`~.Interval`
            The dna sequence interval. Anchoring will expand the
            interval when deletion variants are encountered.

        Returns
        -------
        :class:`list` of :class:`~.VCFVariant`
            A list of variants that possibly affect the given `interval`'s DNA
            sequence.
        """
        if self._mask_indices is not None:
            raise NotImplementedError("masked variants are not yet supported.")
        return _cxx.VCFTable.sequence_variations(self, interval)

    @mock
    def close(self):  # pragma: no cover
        """Close the file handle.

        Note that if you use a `with` statement the file handle is automatically closed::

            # Open the file
            with VCFTable('foo.vcfbin') as vcfbin:
                ...

            # <-- File is now closed

        Notes
        -----
        Closing this table will invalidate all VCFVariants previously accessed via __getitem__.
        """
        mock_unreachable()

    @mock
    @property
    def filename(self):  # pragma: no cover
        """The path to the file from which FORMAT field matrices are retrieved.

        Returns
        -------
        :py:class:`str`
           The path to the VCF file.
        """
        return mock_result(str)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __repr__(self):
        return "<VCFTable, len() = {}>".format(len(self))


########################################################################
