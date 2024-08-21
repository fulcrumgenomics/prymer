"""
# Variant Lookup Class and Methods

This module contains the abstract class
[`VariantLookup`][prymer.api.variant_lookup.VariantLookup] to facilitate retrieval of
variants that overlap a given genomic coordinate range.  Concrete implementations must implement the
[`query()`][prymer.api.variant_lookup.VariantLookup.query] method for retrieving variants
that overlap the given range.

Two concrete implementations are provided that both take a list of VCF files to be queried:

- [`FileBasedVariantLookup`][prymer.api.variant_lookup.FileBasedVariantLookup] -- performs
disk-based retrieval of variants (using a VCF index).  This class is recommended for large VCFs. The
[`disk_based()`][prymer.api.variant_lookup.disk_based] alternative constructor is
provided for easy construction of this object.
- [`VariantOverlapDetector`][prymer.api.variant_lookup.VariantOverlapDetector] -- reads in
variants into memory and uses an
[`pybedlite.overlap_detector.OverlapDetector`](https://pybedlite.readthedocs.io/en/latest/api.html#pybedlite.overlap_detector.OverlapDetector)
for querying.  This class is recommended for small VCFs. The
[`cached()`][prymer.api.variant_lookup.cached] alternative constructor is provided for
easy construction of this object.

Each class can also use minor allele frequency (MAF) to filter variants.

The helper class `SimpleVariant` is included to facilitate VCF querying and reporting out results.

## Examples

```python
>>> from pathlib import Path
>>> lookup = cached(vcf_paths=[Path("./tests/api/data/miniref.variants.vcf.gz")], min_maf=0.00, include_missing_mafs=True)
>>> lookup.query(refname="chr2", start=7999, end=8000)
[SimpleVariant(id='complex-variant-sv-1/1', refname='chr2', pos=8000, ref='T', alt='<DEL>', end=8000, variant_type=<VariantType.OTHER: 'OTHER'>, maf=None)]
>>> variants = lookup.query(refname="chr2", start=7999, end=9900)
>>> len(variants)
14
>>> for v in variants: \
    print(v)
complex-variant-sv-1/1@chr2:8000[T/<DEL> NA]
rare-dbsnp-snp1-1/1@chr2:9000[A/C 0.0010]
common-dbsnp-snp1-1/1@chr2:9010[C/T 0.0100]
common-dbsnp-snp2-1/1@chr2:9020[A/C 0.0200]
rare-ac-an-snp-1/1@chr2:9030[G/A 0.0010]
rare-af-snp-1/1@chr2:9040[C/A 0.0005]
common-ac-an-snp-1/1@chr2:9050[C/G 0.0470]
common-af-snp-1/1@chr2:9060[T/C 0.0400]
common-multiallelic-1/2@chr2:9070[C/A 0.0525]
common-multiallelic-2/2@chr2:9070[C/T 0.0525]
common-insertion-1/1@chr2:9080[A/ACGT 0.0400]
common-deletion-1/1@chr2:9090[CTA/C 0.0400]
common-mixed-1/2@chr2:9100[CA/GG 0.1200]
common-mixed-2/2@chr2:9101[A/ACACA 0.1200]
>>> lookup.query(refname="ch12", start=7999, end=9900)
[]

```
"""  # noqa: E501

import logging
from abc import ABC
from abc import abstractmethod
from dataclasses import dataclass
from dataclasses import field
from enum import auto
from enum import unique
from pathlib import Path
from typing import Optional
from typing import final

import pysam
from fgpyo.vcf import reader
from pybedlite.overlap_detector import Interval
from pybedlite.overlap_detector import OverlapDetector
from pysam import VariantFile
from pysam import VariantRecord
from strenum import UppercaseStrEnum

from prymer.api.span import Span
from prymer.api.span import Strand


@unique
class VariantType(UppercaseStrEnum):
    """Represents the type of variant."""

    SNP = auto()
    MNV = auto()
    INSERTION = auto()
    DELETION = auto()
    OTHER = auto()

    @staticmethod
    def from_alleles(ref: str, alt: str) -> "VariantType":
        """Builds a variant type from the given reference and alternate allele.

        Args:
            ref: the reference allele
            alt: the alternate allele

        Returns:
            the variant type
        """
        variant_type: VariantType
        if "<" in alt:
            variant_type = VariantType.OTHER
        elif (len(ref) == 1) and (len(alt) == 1):
            variant_type = VariantType.SNP
        elif len(ref) == len(alt):
            variant_type = VariantType.MNV
        elif len(ref) < len(alt):
            variant_type = VariantType.INSERTION
        elif len(ref) > len(alt):
            variant_type = VariantType.DELETION
        else:
            raise ValueError(f"Could not determine variant type for ref `{ref}` and alt `{alt}`")
        return variant_type


@dataclass(slots=True, frozen=True, init=True)
class SimpleVariant:
    """Represents a variant from a given genomic range.

    The reference and alternate alleles typically match the VCF from which the variant originated,
    differing when the reference and alternate alleles share common leading bases.  For example,
    given a reference allele of `CACACA` and alternate allele `CACA`, the reference and alternate
    alleles have their common leading bases removed except for a single "anchor base" (i.e. `CAC` is
    removed), yielding the reference allele `ACA` and alternate allele `A`.  This will also change
    the reference variant position (`pos` property), which is defined as the position of the first
    base in the reference allele.

    The `variant_type` is derived from the reference and alternate alleles after the above (see
    [`VariantType.from_alleles()`][prymer.api.variant_lookup.VariantType.from_alleles]).

    Furthermore, the variant end position (`end` property) is defined for insertions as the base
    after the insertion.  For all other variant types, this is the last base in the reference
    allele.

    Attributes:
        id: the variant identifier
        refname: the reference sequence name
        pos: the start position of the variant (1-based inclusive)
        end: the end position of the variant (1-based inclusive)
        ref: the reference base
        alt: the alternate base
        variant_type: the variant type
        maf: optionally, the minor allele frequency
    """

    id: str
    refname: str
    pos: int
    ref: str
    alt: str
    end: int = field(init=False)
    variant_type: VariantType = field(init=False)
    maf: Optional[float] = None

    def __post_init__(self) -> None:
        # Simplify the ref/alt alleles
        # find the leading prefix of the ref and alternate allele
        prefix_length = 0
        for r, a in zip(self.ref, self.alt, strict=False):
            if r != a:
                break
            prefix_length += 1
        if prefix_length > 1:  # > 1 to always keep an anchor base
            offset = prefix_length - 1
            object.__setattr__(self, "ref", self.ref[offset:])
            object.__setattr__(self, "alt", self.alt[offset:])
            object.__setattr__(self, "pos", self.pos + offset)

        # Derive the end and variant_type from the ref and alt alleles.
        object.__setattr__(
            self, "variant_type", VariantType.from_alleles(ref=self.ref, alt=self.alt)
        )

        # span all reference bases (e.g. SNP -> 1bp, MNV -> len(MNV), DEL -> len(DEL),
        end: int = self.pos + len(self.ref) - 1
        if self.variant_type == VariantType.INSERTION:
            # span the base preceding and following the insertion, so add a base at the end
            end += 1
        object.__setattr__(self, "end", end)

    def __str__(self) -> str:
        """Compact String representation of the variant that includes all relevant info."""
        maf_string = f"{self.maf:.4f}" if self.maf is not None else "NA"
        return f"{self.id}@{self.refname}:{self.pos}[{self.ref}/{self.alt} {maf_string}]"

    def to_span(self) -> Span:
        """Creates a Span object that represents the genomic span of this variant.
        Insertions will span the base preceding and following the inserted bases."""
        return Span(refname=self.refname, start=self.pos, end=self.end, strand=Strand.POSITIVE)

    @staticmethod
    def build(variant: VariantRecord) -> list["SimpleVariant"]:
        """Convert `pysam.VariantRecord` to `SimpleVariant`. Only the first ALT allele is used."""
        maf = calc_maf_from_filter(variant)
        simple_variants: list[SimpleVariant] = []

        for i, alt in enumerate(variant.alts, start=1):
            simple_variant = SimpleVariant(
                id=f"{variant.id}-{i}/{len(variant.alts)}",
                refname=variant.chrom,
                pos=variant.pos,
                ref=variant.ref,
                alt=alt,
                maf=maf,
            )
            simple_variants.append(simple_variant)

        return simple_variants


@dataclass(slots=True, frozen=True, init=True)
class _VariantInterval(Interval):
    """Intermediate class that facilitates use of pybedlite.overlap_detector with a `SimpleVariant`.
    This must be an @attr.s because the `pybedlite.Interval` is an @attr.s."""

    refname: str
    start: int
    end: int
    variant: SimpleVariant
    negative: bool = False
    name: Optional[str] = None

    @staticmethod
    def build(simple_variant: SimpleVariant) -> "_VariantInterval":
        return _VariantInterval(
            refname=simple_variant.refname,
            start=simple_variant.pos - 1,  # pybedlite is 0-based open ended
            end=simple_variant.pos + len(simple_variant.ref) - 1,
            variant=simple_variant,
        )


class VariantLookup(ABC):
    """Base class to represent a variant from a given genomic range.

    Attributes:
        vcf_paths: the paths to the source VCFs for the variants
        min_maf: optionally, return only variants with at least this minor allele frequency
        include_missing_mafs: when filtering variants with a minor allele frequency,
            `True` to include variants with no annotated minor allele frequency, otherwise `False`.
            If no minor allele frequency is given, then this parameter does nothing.
    """

    def __init__(
        self,
        vcf_paths: list[Path],
        min_maf: Optional[float],
        include_missing_mafs: bool,
    ) -> None:
        self.vcf_paths: list[Path] = vcf_paths
        self.min_maf: Optional[float] = min_maf
        self.include_missing_mafs: bool = include_missing_mafs

    @final
    def query(
        self,
        refname: str,
        start: int,
        end: int,
        maf: Optional[float] = None,
        include_missing_mafs: bool = None,
    ) -> list[SimpleVariant]:
        """Gets all variants that overlap a genomic range, optionally filter based on MAF threshold.

        Args:
            refname: the reference name
            start: the 1-based start position
            end: the 1-based end position
            maf: the MAF of the variant
            include_missing_mafs: whether to include variants with a missing MAF
                (overrides self.include_missing_mafs)
        """
        if maf is None:
            maf = self.min_maf
        if include_missing_mafs is None:
            include_missing_mafs = self.include_missing_mafs

        variants = self._query(refname=refname, start=start, end=end)
        if len(variants) == 0:
            logging.debug(f"No variants extracted from region of interest: {refname}:{start}-{end}")
        if maf is None or maf <= 0.0:
            return variants
        elif include_missing_mafs:  # return variants with a MAF above threshold or missing
            return [v for v in variants if (v.maf is None or v.maf >= maf)]
        else:
            return [v for v in variants if v.maf is not None and v.maf >= maf]

    @staticmethod
    def to_variants(
        variants: list[VariantRecord], source_vcf: Optional[Path] = None
    ) -> list[SimpleVariant]:
        """Converts a list of `pysam.VariantRecords` to a list of `SimpleVariants` for ease of use.
        Filters variants based on their FILTER status, and sorts by start position.

        Args:
            variants: the variants to convert
            source_vcf: the optional path to the source VCF, used only for debug messages

        Returns:
            a list of `SimpleVariants`, one per alternate allele per variant.
        """
        simple_vars = []
        for variant in variants:
            if (
                "PASS" in list(variant.filter) or len(list(variant.filter)) == 0
            ):  # if passing or empty filters
                simple_variants = SimpleVariant.build(variant)
                if any(v.variant_type == VariantType.OTHER for v in simple_variants):
                    logging.debug(
                        f"Input VCF file {source_vcf} may contain complex variants: {variant}"
                    )
                simple_vars.extend(simple_variants)
        return sorted(simple_vars, key=lambda v: (v.pos, v.id))

    @abstractmethod
    def _query(self, refname: str, start: int, end: int) -> list[SimpleVariant]:
        """Subclasses must implement this method."""


class FileBasedVariantLookup(VariantLookup):
    """Implementation of VariantLookup that queries against indexed VCF files each time a query is
    performed. Assumes the index is located adjacent to the VCF file and has the same base name with
     either a .csi or .tbi suffix."""

    def __init__(self, vcf_paths: list[Path], min_maf: Optional[float], include_missing_mafs: bool):
        self._readers: list[VariantFile] = []
        super().__init__(
            vcf_paths=vcf_paths, min_maf=min_maf, include_missing_mafs=include_missing_mafs
        )
        if len(vcf_paths) == 0:
            raise ValueError("No VCF paths given to query.")
        for path in vcf_paths:
            if (
                not path.with_suffix(path.suffix + ".csi").is_file()
                and not path.with_suffix(path.suffix + ".tbi").is_file()
            ):
                raise ValueError(f"Cannot perform fetch with missing index file for VCF: {path}")
            open_fh = pysam.VariantFile(str(path))
            self._readers.append(open_fh)

    def _query(self, refname: str, start: int, end: int) -> list[SimpleVariant]:
        """Queries variants from the VCFs used by this lookup and returns a `SimpleVariant`."""
        simple_variants: list[SimpleVariant] = []
        for fh, path in zip(self._readers, self.vcf_paths, strict=True):
            if not fh.header.contigs.get(refname):
                logging.debug(f"Header in VCF file {path} does not contain chromosome {refname}.")
                continue
            #  pysam.fetch is 0-based, half-open
            variants = [variant for variant in fh.fetch(contig=refname, start=start - 1, end=end)]
            simple_variants.extend(self.to_variants(variants, source_vcf=path))
        return sorted(simple_variants, key=lambda x: x.pos)


class VariantOverlapDetector(VariantLookup):
    """Implements `VariantLookup` by reading the entire VCF into memory and loading the resulting
    Variants into an `OverlapDetector`."""

    def __init__(self, vcf_paths: list[Path], min_maf: Optional[float], include_missing_mafs: bool):
        super().__init__(
            vcf_paths=vcf_paths, min_maf=min_maf, include_missing_mafs=include_missing_mafs
        )
        self._overlap_detector: OverlapDetector[_VariantInterval] = OverlapDetector()

        if len(vcf_paths) == 0:
            raise ValueError("No VCF paths given to query.")

        for path in vcf_paths:
            if (
                not path.with_suffix(path.suffix + ".csi").is_file()
                and not path.with_suffix(path.suffix + ".tbi").is_file()
            ):
                raise ValueError(f"Cannot perform fetch with missing index file for VCF: {path}")

            with reader(path) as fh:
                variant_intervals = iter(
                    _VariantInterval.build(simple_variant)
                    for variant in fh
                    for simple_variant in self.to_variants([variant], source_vcf=path)
                )
                self._overlap_detector.add_all(variant_intervals)

    def _query(self, refname: str, start: int, end: int) -> list[SimpleVariant]:
        """Queries variants from the VCFs used by this lookup."""
        query = Interval(
            refname=refname, start=start - 1, end=end
        )  # interval is half-open, 0-based
        overlapping_variants = (
            variant_interval.variant
            for variant_interval in self._overlap_detector.get_overlaps(query)
        )
        return sorted(overlapping_variants, key=lambda v: (v.pos, v.id))


# module-level functions
def calc_maf_from_filter(variant: pysam.VariantRecord) -> Optional[float]:
    """Calculates minor allele frequency (MAF) based on whether `VariantRecord` denotes
    CAF, AF, AC and AN in the INFO field. In none of those are found, looks in the FORMAT
    field for genotypes and calculates MAF.

    If the variant has multiple _alternate_ alleles, then the returned MAF is the sum of MAFs across
    all alternate alleles.

    Args:
        variant: the variant from which the MAF will be computed

    Returns:
        the minor allele frequency (MAF)
    """

    maf: Optional[float] = None
    if "CAF" in variant.info:
        # Assumes CAF is a list of allele frequencies, with the first being the reference allele
        maf = 1 - float(variant.info["CAF"][0])
    elif "AF" in variant.info:
        maf = sum(float(af) for af in variant.info["AF"])
    elif "AC" in variant.info and "AN" in variant.info:
        ac = sum(int(ac) for ac in variant.info["AC"])
        an = int(variant.info["AN"])
        maf = ac / an
    elif len(list(variant.samples)) > 0:  # if genotypes are not empty
        gts = [idx for sample in variant.samples.values() for idx in sample["GT"] if "GT" in sample]
        if len(gts) > 0:
            num_alt = sum(1 for idx in gts if idx != 0)
            maf = num_alt / len(gts)

    return maf


def cached(
    vcf_paths: list[Path], min_maf: float, include_missing_mafs: bool = False
) -> VariantOverlapDetector:
    """Constructs a `VariantLookup` that caches all variants in memory for fast lookup.
    Appropriate for small VCFs."""
    return VariantOverlapDetector(
        vcf_paths=vcf_paths, min_maf=min_maf, include_missing_mafs=include_missing_mafs
    )


def disk_based(
    vcf_paths: list[Path], min_maf: float, include_missing_mafs: bool = False
) -> FileBasedVariantLookup:
    """Constructs a `VariantLookup` that queries indexed VCFs on disk for each lookup.
    Appropriate for large VCFs."""
    return FileBasedVariantLookup(
        vcf_paths=vcf_paths, min_maf=min_maf, include_missing_mafs=include_missing_mafs
    )
