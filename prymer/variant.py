import logging
from dataclasses import dataclass
from pathlib import Path
from types import TracebackType
from typing import ContextManager
from typing import Optional

import pysam
from pysam import VariantFile
from pysam import VariantRecord

from prymer.model import Span
from prymer.model import Strand

_logger = logging.getLogger(__name__)


@dataclass(slots=True, frozen=True, init=True)
class SimpleVariant:
    """Represents a variant from a given genomic range.

    Each variant has a single reference and alt allele.  Start and end positions are both 1-based
    inclusive.  Unlike VCF, empty alleles are allowed and positions are given as follows:

    - SNPs and MNPS will have start/end set to the first and last substituted bases
    - Deletions will have start/end set to the first and last _deleted_ base
    - Insertions will have start/end set to the bases immediately before and after the insertion

    Attributes:
        id: the variant identifier
        refname: the reference sequence name
        start: the start position of the variant (1-based inclusive)
        end: the end position of the variant (1-based inclusive)
        ref: the reference base
        alt: the alternate base
        af: the alt allele frequency if available
    """

    id: str
    refname: str
    start: int
    end: int
    ref: str
    alt: str
    af: Optional[float] = None

    def __str__(self) -> str:
        """Compact String representation of the variant that includes all relevant info."""
        af_string = f"{self.af:.4f}" if self.af is not None else "NA"
        return f"{self.id}@{self.refname}:{self.start}[{self.ref}/{self.alt} {af_string}]"

    def to_span(self) -> Span:
        """Creates a Span object that represents the genomic span of this variant.
        Insertions will span the base preceding and following the inserted bases."""
        return Span(refname=self.refname, start=self.start, end=self.end, strand=Strand.POSITIVE)

    @staticmethod
    def build(variant: VariantRecord) -> list["SimpleVariant"]:
        """
        Convert `pysam.VariantRecord` to one or more `SimpleVariant`s.
        """
        # Skip symbolic ref alleles
        if "<" in variant.ref:
            return []

        ref = variant.ref
        alt_afs = SimpleVariant._extract_mafs(variant)
        if alt_afs is not None and len(alt_afs) != len(variant.alts):
            raise ValueError(f"Different number of allele frequencies vs. alt alleles: {variant}")

        simple_variants: list[SimpleVariant] = []

        for i, alt in enumerate(variant.alts):
            # Skip symbolic alt alleles
            if "<" in alt:
                continue

            if len(alt) < len(ref) and ref.startswith(alt):
                # Deletion
                padding_bases = len(alt)
                deleted_bases = len(ref) - padding_bases
                start = variant.pos + padding_bases
                end = start + deleted_bases - 1
                ref = ref[padding_bases:]
                alt = ""
            elif len(alt) > len(ref) and alt.startswith(ref):
                # Insertion
                padding_bases = len(ref)
                deleted_bases = len(ref) - padding_bases
                start = variant.pos + padding_bases - 1
                end = start + 1
                ref = ""
                alt = alt[padding_bases:]
            else:
                # Substitution (simple or complex)
                start = variant.pos
                end = start + len(ref) - 1

            id = variant.id or f"{variant.chrom}:{start}"
            if len(variant.alts) > 1:
                id = f"{id}-{i + 1}/{len(variant.alts)}"

            simple_variant = SimpleVariant(
                id=id,
                refname=variant.chrom,
                start=start,
                end=end,
                ref=ref,
                alt=alt,
                af=alt_afs[i] if alt_afs is not None else None,
            )
            simple_variants.append(simple_variant)

        return simple_variants

    @staticmethod
    def _extract_mafs(variant: pysam.VariantRecord) -> Optional[list[float]]:
        """
        Calculates/extracts alt allele frequencies (MAF) based on whether `VariantRecord` has
        CAF, AF, AC and AN in the INFO field.

        Args:
            variant: the variant from which the MAF will be computed

        Returns:
            the list of frequencies of alternate alleles
        """
        if "CAF" in variant.info:
            # Assumes CAF is a list of allele frequencies, with the first being the reference allele
            return [float(af) for af in variant.info["CAF"][1:]]
        elif "AF" in variant.info:
            return [float(af) for af in variant.info["AF"]]
        elif "AC" in variant.info and "AN" in variant.info:
            an = float(variant.info["AN"])
            return [float(ac) / an for ac in variant.info["AC"]]

        return None


class VariantLookup(ContextManager):
    """Class to provide lookup and filtering of variants from one or more VCF files. Additionally
    provides methods for masking fragments of sequence for variant bases.

    Attributes:
        vcf_paths: the paths to the source VCFs for the variants
        min_maf: only use variants with an allele frequency above the minimum given.
        include_missing_mafs: when filtering variants by allele frequency should variants whose
          allele frequency cannot be determined be included?  Defaults to true when `min_maf` is 0
          otherwise false.
    """

    def __init__(
        self, vcf_paths: list[Path], min_maf: float = 0, include_missing_mafs: Optional[bool] = None
    ) -> None:
        if len(vcf_paths) == 0:
            raise ValueError("No VCF paths given to query.")

        self.vcf_paths: list[Path] = vcf_paths
        self.min_maf: float = max(0.0, min_maf)
        self.include_missing_mafs: bool = (
            include_missing_mafs if include_missing_mafs is not None else min_maf == 0.0
        )
        self._readers: list[VariantFile] = []

        for path in vcf_paths:
            if (
                not path.with_suffix(path.suffix + ".csi").is_file()
                and not path.with_suffix(path.suffix + ".tbi").is_file()
            ):
                raise ValueError(f"Cannot find index file for VCF: {path}")

            open_fh = pysam.VariantFile(str(path))
            self._readers.append(open_fh)

    def __enter__(self) -> "VariantLookup":
        """Enter the context manager."""
        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        """Exit this context manager while closing the underlying VCF handles."""
        self.close()
        return None

    def query(self, refname: str, start: int, end: int) -> list[SimpleVariant]:
        """Fetches variants that overlap a genomic range.

        Args:
            refname: the reference name
            start: the 1-based start position
            end: the 1-based inclusive end position
        """
        variants: list[SimpleVariant] = []

        for fh, path in zip(self._readers, self.vcf_paths, strict=True):
            if not fh.header.contigs.get(refname):
                _logger.debug(f"Header in VCF file {path} does not contain chromosome {refname}.")
                continue

            #  pysam.fetch is 0-based, half-open
            for rec in fh.fetch(contig=refname, start=start - 1, end=end):
                # Skip over variants with a non-PASS filter
                filters = list(rec.filter)
                if len(filters) > 0 and "PASS" not in filters:
                    continue

                for variant in SimpleVariant.build(rec):
                    if variant.af is None and self.include_missing_mafs:
                        variants.append(variant)
                    elif variant.af is not None and variant.af >= self.min_maf:
                        variants.append(variant)

        variants.sort(key=lambda v: v.start)
        return variants

    def mask(self, bases: str, region: Span) -> tuple[str, str]:
        """
        Takes a sequence and returns soft and hard masked version of it.  Bases in the sequence
        are masked if a variant overlaps that base and passes frequency filtering.

        The input sequence is upper-cased, and then masked bases are converted to lower case in the
        soft-masked version and to Ns in the hard-masked version.

        Args:
            bases: the bases to be masked
            region: a span describing the genomic region from which the bases are drawn

        Returns:
            A tuple of two sequences where the first is soft-masked and the second is hard masked
        """
        if len(bases) != region.length:
            raise ValueError("Region and bases have different lengths: {region} vs. {bases}")

        bases = bases.upper()
        variants = self.query(refname=region.refname, start=region.start, end=region.end)

        if len(variants) == 0:
            return bases, bases
        else:
            soft = list(bases)
            hard = list(bases)

            for v in variants:
                for pos in range(max(v.start, region.start), min(v.end, region.end) + 1):
                    idx = pos - region.start
                    soft[idx] = soft[idx].lower()
                    hard[idx] = "N"

            return "".join(soft), "".join(hard)

    def close(self) -> None:
        """Close the underlying VCF file handles."""
        for handle in self._readers:
            handle.close()
