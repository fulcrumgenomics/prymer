import logging
import random
from contextlib import closing
from dataclasses import dataclass
from dataclasses import replace
from pathlib import Path
from typing import Optional

import fgpyo.vcf.builder
import pytest
from fgpyo.vcf.builder import VariantBuilder
from fgpyo.vcf.builder import VcfFieldNumber
from fgpyo.vcf.builder import VcfFieldType
from pysam import VariantRecord

from prymer.api.span import Span
from prymer.api.span import Strand
from prymer.api.variant_lookup import FileBasedVariantLookup
from prymer.api.variant_lookup import SimpleVariant
from prymer.api.variant_lookup import VariantOverlapDetector
from prymer.api.variant_lookup import VariantType
from prymer.api.variant_lookup import cached
from prymer.api.variant_lookup import calc_maf_from_filter
from prymer.api.variant_lookup import disk_based


@pytest.mark.parametrize(
    "ref, alt, variant_type",
    [
        ("A", "C", VariantType.SNP),
        ("AA", "CC", VariantType.MNV),
        ("ACAC", "ACAC", VariantType.MNV),
        ("A", "ACAC", VariantType.INSERTION),
        ("AC", "ACAC", VariantType.INSERTION),
        ("ACA", "ACAC", VariantType.INSERTION),
        ("ACAC", "ACA", VariantType.DELETION),
        ("ACAC", "AC", VariantType.DELETION),
        ("ACAC", "A", VariantType.DELETION),
        ("A", "<DEL>", VariantType.OTHER),
    ],
)
def test_variant_type_build(ref: str, alt: str, variant_type: VariantType) -> None:
    assert variant_type == VariantType.from_alleles(ref=ref, alt=alt)


@pytest.fixture(scope="session")
def vcf_path(data_dir: Path) -> Path:
    """Test VCF data: returns a path to a VCF containing 12 variants.
    8 SNPs, 1 deletion, 1 insertion, 1 mixed variant, and 1 complex variant."""
    return Path(__file__).parent / "data" / "miniref.variants.vcf.gz"


@pytest.fixture
def temp_missing_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Invalid path to test VCF."""
    return tmp_path_factory.mktemp("test_missing_vcf")


@pytest.fixture()
def mini_chr1_vcf() -> Path:
    """Create an in-memory mini VCF using `fgpyo.variant_builder` for testing. These variants
    are different from what is contained in the miniref.variants.vcf.gz file."""
    return Path(__file__).parent / "data" / "mini_chr1.vcf.gz"


@pytest.fixture()
def mini_chr3_vcf() -> Path:
    """Create an in-memory mini VCF using `fgpyo.variant_builder` for testing. These variants
    are different from what is contained in the miniref.variants.vcf.gz file."""
    return Path(__file__).parent / "data" / "mini_chr3.vcf.gz"


@pytest.fixture()
def sample_vcf() -> list[VariantRecord]:
    """Create an in-memory VCF using `fgpyo.variant_builder` for testing.
    These variants here are the same as what is contained in the miniref.variants.vcf.gz file."""
    variant_builder = VariantBuilder(
        sd={"chr1": {"ID": "chr1", "length": 577}, "chr2": {"ID": "chr2", "length": 9821}}
    )
    variant_builder.add_info_header(name="AC", field_type=VcfFieldType.INTEGER, number=2)
    variant_builder.add_info_header(
        name="AF", field_type=VcfFieldType.FLOAT, number=VcfFieldNumber.NUM_ALT_ALLELES
    )
    variant_builder.add_info_header(name="AN", field_type=VcfFieldType.INTEGER, number=1)
    variant_builder.add_info_header(name="CAF", field_type=VcfFieldType.STRING, number=1)
    variant_builder.add_info_header(name="COMMON", field_type=VcfFieldType.INTEGER)
    variant_builder.add_info_header(name="SVTYPE", field_type=VcfFieldType.STRING)
    variant_builder.add_info_header(name="SVLEN", field_type=VcfFieldType.INTEGER)
    variant_builder.add(
        contig="chr2",
        pos=8000,
        id="complex-variant-sv",
        ref="T",
        alts="<DEL>",
        info={"SVTYPE": "DEL", "SVLEN": -105},
    )
    variant_builder.add(
        contig="chr2",
        pos=9000,
        id="rare-dbsnp-snp1",
        ref="A",
        alts="C",
        info={"CAF": "0.999,0.001", "COMMON": 0},
    )
    variant_builder.add(
        contig="chr2",
        pos=9010,
        id="common-dbsnp-snp1",
        ref="C",
        alts="T",
        info={"CAF": "0.99,0.01", "COMMON": 1},
    )
    variant_builder.add(
        contig="chr2",
        pos=9020,
        id="common-dbsnp-snp2",
        ref="A",
        alts="C",
        info={"CAF": "0.98,0.02", "COMMON": 0},
    )
    variant_builder.add(
        contig="chr2",
        pos=9030,
        id="rare-ac-an-snp",
        ref="G",
        alts="A",
        info={"AC": [1, 0], "AN": 1000},
    )
    variant_builder.add(
        contig="chr2", pos=9040, id="rare-af-snp", ref="C", alts="A", info={"AF": 0.0004815}
    )
    variant_builder.add(
        contig="chr2",
        pos=9050,
        id="common-ac-an-snp",
        ref="C",
        alts="G",
        info={"AC": [47, 0], "AN": 1000},
    )
    variant_builder.add(
        contig="chr2", pos=9060, id="common-af-snp", ref="T", alts="C", info={"AF": 0.04}
    )
    variant_builder.add(
        contig="chr2",
        pos=9070,
        id="common-multiallelic",
        ref="C",
        alts=("A", "T"),
        info={"AC": [3, 18], "AN": 400},
    )
    variant_builder.add(
        contig="chr2", pos=9080, id="common-insertion", ref="A", alts="ACGT", info={"AF": 0.04}
    )
    variant_builder.add(
        contig="chr2", pos=9090, id="common-deletion", ref="CTA", alts="C", info={"AF": 0.04}
    )
    variant_builder.add(
        contig="chr2",
        pos=9100,
        id="common-mixed",
        ref="CA",
        alts=("GG", "CACACA"),
        info={"AF": [0.04, 0.08]},
    )

    return variant_builder.to_sorted_list()


@dataclass(init=True, frozen=True)
class SimpleVariantTestCase:
    """Test case for a `SimpleVariant`.

    Attributes:
        simple_variant: the `SimpleVariant` to test
        str_rep: the expected string representation
        span: the expected span
        var_type: the expected variant type
    """

    simple_variant: SimpleVariant
    str_rep: str
    span: Span
    var_type: VariantType


def build_simple_variant_test_cases() -> list[SimpleVariantTestCase]:
    """Builds a list of test cases for `SimpleVariant` methods."""
    return [
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="complex-variant-sv-1/1",
                refname="chr2",
                pos=8000,
                ref="T",
                alt="<DEL>",
                maf=None,
            ),
            str_rep="complex-variant-sv-1/1@chr2:8000[T/<DEL> NA]",
            span=Span(refname="chr2", start=8000, end=8000, strand=Strand.POSITIVE),
            var_type=VariantType.OTHER,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="rare-dbsnp-snp1-1/1",
                refname="chr2",
                pos=9000,
                ref="A",
                alt="C",
                maf=0.001,
            ),
            str_rep="rare-dbsnp-snp1-1/1@chr2:9000[A/C 0.0010]",
            span=Span(refname="chr2", start=9000, end=9000, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            SimpleVariant(
                id="common-dbsnp-snp1-1/1",
                refname="chr2",
                pos=9010,
                ref="C",
                alt="T",
                maf=0.01,
            ),
            str_rep="common-dbsnp-snp1-1/1@chr2:9010[C/T 0.0100]",
            span=Span(refname="chr2", start=9010, end=9010, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-dbsnp-snp2-1/1",
                refname="chr2",
                pos=9020,
                ref="A",
                alt="C",
                maf=0.02,
            ),
            str_rep="common-dbsnp-snp2-1/1@chr2:9020[A/C 0.0200]",
            span=Span(refname="chr2", start=9020, end=9020, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="rare-ac-an-snp-1/1", refname="chr2", pos=9030, ref="G", alt="A", maf=0.001
            ),
            str_rep="rare-ac-an-snp-1/1@chr2:9030[G/A 0.0010]",
            span=Span(refname="chr2", start=9030, end=9030, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="rare-af-snp-1/1",
                refname="chr2",
                pos=9040,
                ref="C",
                alt="A",
                maf=0.00048149999929592013,
            ),
            str_rep="rare-af-snp-1/1@chr2:9040[C/A 0.0005]",
            span=Span(refname="chr2", start=9040, end=9040, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-ac-an-snp-1/1", refname="chr2", pos=9050, ref="C", alt="G", maf=0.047
            ),
            str_rep="common-ac-an-snp-1/1@chr2:9050[C/G 0.0470]",
            span=Span(refname="chr2", start=9050, end=9050, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-af-snp-1/1",
                refname="chr2",
                pos=9060,
                ref="T",
                alt="C",
                maf=0.03999999910593033,
            ),
            str_rep="common-af-snp-1/1@chr2:9060[T/C 0.0400]",
            span=Span(refname="chr2", start=9060, end=9060, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-multiallelic-1/2", refname="chr2", pos=9070, ref="C", alt="A", maf=0.0525
            ),
            str_rep="common-multiallelic-1/2@chr2:9070[C/A 0.0525]",
            span=Span(refname="chr2", start=9070, end=9070, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-multiallelic-2/2", refname="chr2", pos=9070, ref="C", alt="T", maf=0.0525
            ),
            str_rep="common-multiallelic-2/2@chr2:9070[C/T 0.0525]",
            span=Span(refname="chr2", start=9070, end=9070, strand=Strand.POSITIVE),
            var_type=VariantType.SNP,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-insertion-1/1",
                refname="chr2",
                pos=9080,
                ref="A",
                alt="ACGT",
                maf=0.03999999910593033,
            ),
            str_rep="common-insertion-1/1@chr2:9080[A/ACGT 0.0400]",
            span=Span(refname="chr2", start=9080, end=9081, strand=Strand.POSITIVE),
            var_type=VariantType.INSERTION,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-deletion-1/1",
                refname="chr2",
                pos=9090,
                ref="CTA",
                alt="C",
                maf=0.03999999910593033,
            ),
            str_rep="common-deletion-1/1@chr2:9090[CTA/C 0.0400]",
            span=Span(refname="chr2", start=9090, end=9092, strand=Strand.POSITIVE),  # end adjusted
            var_type=VariantType.DELETION,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-mixed-1/2",
                refname="chr2",
                pos=9100,
                ref="CA",
                alt="GG",
                maf=0.12,
            ),
            str_rep="common-mixed-1/2@chr2:9100[CA/GG 0.1200]",
            span=Span(refname="chr2", start=9100, end=9101, strand=Strand.POSITIVE),
            var_type=VariantType.MNV,
        ),
        SimpleVariantTestCase(
            simple_variant=SimpleVariant(
                id="common-mixed-2/2",
                refname="chr2",
                pos=9101,
                ref="A",
                alt="ACACA",
                maf=0.12,
            ),
            str_rep="common-mixed-2/2@chr2:9101[A/ACACA 0.1200]",
            span=Span(refname="chr2", start=9101, end=9102, strand=Strand.POSITIVE),
            var_type=VariantType.INSERTION,
        ),
    ]


def _round_simple_variant(simple_variant: SimpleVariant, digits: int = 5) -> SimpleVariant:
    """Returns a copy of the simple variant after rounding the minor allele frequency
    (maf attribute) to five decimals"""
    maf: Optional[float] = None
    if simple_variant.maf is not None:
        maf = round(simple_variant.maf, digits)
    return replace(simple_variant, maf=maf)


VALID_SIMPLE_VARIANT_TEST_CASES: list[SimpleVariantTestCase] = build_simple_variant_test_cases()
"""Test cases for `SimpleVariant` corresponding to the variants in the fixture `vcf_path`"""

VALID_SIMPLE_VARIANTS_APPROX: list[SimpleVariant] = [
    _round_simple_variant(test_case.simple_variant) for test_case in VALID_SIMPLE_VARIANT_TEST_CASES
]
"""`SimpleVariant`s corresponding to VALID_SIMPLE_VARIANT_TEST_CASES but rounded to five digits"""


def get_simple_variant_approx_by_id(*variant_id: str) -> list[SimpleVariant]:
    return [
        next(
            simple_variant
            for simple_variant in VALID_SIMPLE_VARIANTS_APPROX
            if simple_variant.id == _id
        )
        for _id in variant_id
    ]


def variant_overlap_detector_query(
    detector: VariantOverlapDetector,
    refname: str,
    start: int,
    end: int,
    maf: Optional[float] = None,
    include_missing_mafs: bool = None,
) -> list[SimpleVariant]:
    return [
        _round_simple_variant(variant)
        for variant in detector.query(
            refname=refname,
            start=start,
            end=end,
            maf=maf,
            include_missing_mafs=include_missing_mafs,
        )
    ]


@pytest.mark.parametrize("test_case", VALID_SIMPLE_VARIANT_TEST_CASES)
def test_simple_variant_var_type(test_case: SimpleVariantTestCase) -> None:
    """Test the `SimpleVariant.variant_type` property."""
    assert test_case.simple_variant.variant_type == test_case.var_type


@pytest.mark.parametrize("test_case", VALID_SIMPLE_VARIANT_TEST_CASES)
def test_simple_variant_str_repr(test_case: SimpleVariantTestCase) -> None:
    """Test the `SimpleVariant.__str__` property."""
    assert test_case.simple_variant.__str__() == test_case.str_rep


@pytest.mark.parametrize("test_case", VALID_SIMPLE_VARIANT_TEST_CASES)
def test_simple_variant_to_span(test_case: SimpleVariantTestCase) -> None:
    """Test the `SimpleVariant.span` property."""
    assert test_case.simple_variant.to_span() == test_case.span


def test_simple_variant_conversion(vcf_path: Path, sample_vcf: list[VariantRecord]) -> None:
    """Test that `pysam.VariantRecords` are converted properly to `SimpleVariants`.
    We pass a valid test data path here to leverage an existing .vcf.gz.tbi file,
    which is required for class instantiation.
    We use the in-memory `VariantBuilder` here to keep test data consistent."""

    variant_overlap_detector = VariantOverlapDetector(
        vcf_paths=[vcf_path], min_maf=0.0, include_missing_mafs=True
    )
    # overcome rounding differences
    actual_simple_variants = [
        _round_simple_variant(v) for v in variant_overlap_detector.to_variants(sample_vcf, vcf_path)
    ]
    assert actual_simple_variants == VALID_SIMPLE_VARIANTS_APPROX


def test_simple_variant_conversion_logs_file_based(
    vcf_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Test that `to_variants()` logs a debug message with no pysam.VariantRecords to convert."""
    caplog.set_level(logging.DEBUG)
    with FileBasedVariantLookup(
        vcf_paths=[vcf_path], min_maf=0.01, include_missing_mafs=False
    ) as variant_lookup:
        variant_lookup.query(refname="foo", start=1, end=2)
        assert "No variants extracted from region of interest" in caplog.text


def test_simple_variant_conversion_logs_non_file_based(
    vcf_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Test that `to_variants()` logs a debug message with no pysam.VariantRecords to convert."""
    caplog.set_level(logging.DEBUG)
    variant_lookup = VariantOverlapDetector(
        vcf_paths=[vcf_path], min_maf=0.01, include_missing_mafs=False
    )
    variant_lookup.query(refname="foo", start=1, end=2)
    assert "No variants extracted from region of interest" in caplog.text


def test_missing_index_file_raises(temp_missing_path: Path) -> None:
    """Test that both VariantLookup objects raise an error with a missing index file."""
    with pytest.raises(ValueError, match="Cannot perform fetch with missing index file for VCF"):
        with closing(
            disk_based(vcf_paths=[temp_missing_path], min_maf=0.01, include_missing_mafs=False)
        ):
            pass
    with pytest.raises(ValueError, match="Cannot perform fetch with missing index file for VCF"):
        cached(vcf_paths=[temp_missing_path], min_maf=0.01, include_missing_mafs=False)


def test_missing_vcf_files_raises() -> None:
    """Test that an error is raised when no VCF_paths are provided."""
    with pytest.raises(ValueError, match="No VCF paths given to query"):
        with closing(disk_based(vcf_paths=[], min_maf=0.01, include_missing_mafs=False)):
            pass
    with pytest.raises(ValueError, match="No VCF paths given to query"):
        cached(vcf_paths=[], min_maf=0.01, include_missing_mafs=False)


@pytest.mark.parametrize("random_seed", [1, 10, 100, 1000, 10000])
def test_vcf_header_missing_chrom(
    mini_chr1_vcf: Path,
    mini_chr3_vcf: Path,
    vcf_path: Path,
    caplog: pytest.LogCaptureFixture,
    random_seed: int,
) -> None:
    """Test whether a missing chromosome of interest from one VCF (in a list of VCFs) will
    log a debug message.
    `vcf_path` contains only variants from chr2
    `mini_chr1_vcf` contains only 3 variants from chr1 (at positions 200, 300, and 400),
    `mini_chr3_vcf` contains only 3 variants from chr3 (at positions 6000, 6010, and 6020)."""
    caplog.set_level(logging.DEBUG)
    vcf_paths = [vcf_path, mini_chr1_vcf, mini_chr3_vcf]
    random.Random(random_seed).shuffle(vcf_paths)
    with FileBasedVariantLookup(
        vcf_paths=vcf_paths, min_maf=0.00, include_missing_mafs=True
    ) as variant_lookup:
        variants_of_interest = variant_lookup.query(
            refname="chr2", start=7999, end=9900
        )  # (chr2 only in vcf_path)
    # Should find all 12 variants from vcf_path (no filtering), with two variants having two
    # alternate alleles
    assert len(variants_of_interest) == 14
    expected_error_msg = "does not contain chromosome"
    assert expected_error_msg in caplog.text


@pytest.mark.parametrize("test_case", VALID_SIMPLE_VARIANT_TEST_CASES)
def test_calc_maf_from_filter(
    test_case: SimpleVariantTestCase, sample_vcf: list[VariantRecord]
) -> None:
    """Test that calculating MAF function is working as expected. Importantly, VALID_SIMPLE_VARIANTS
    are in the same order as the records in the sample_vcf."""
    calculated_mafs: list[float] = []
    expected_mafs = [testcase.simple_variant.maf for testcase in VALID_SIMPLE_VARIANT_TEST_CASES]

    for record in sample_vcf:
        maf = calc_maf_from_filter(record)
        calculated_mafs.extend(maf for _ in record.alts)
    assert calculated_mafs == pytest.approx(expected_mafs)


def test_calc_maf_from_gt_only() -> None:
    """Test that calc_maf using GTs works in the explicit absence of INFO."""
    variant_builder = fgpyo.vcf.builder.VariantBuilder(
        sample_ids=["sample_0", "sample_1", "sample_2"]
    )
    variant_builder._build_header_string()
    variant_builder.add(
        samples={"sample_0": {"GT": (1, 0)}, "sample_1": {"GT": (1, 1)}, "sample_2": {"GT": (0, 0)}}
    )
    for rec in variant_builder.to_sorted_list():
        assert calc_maf_from_filter(rec) == 0.5


def test_variant_overlap_detector_query(vcf_path: Path) -> None:
    """Test `VariantOverlapDetector.query()` positional filtering."""
    variant_overlap_detector = VariantOverlapDetector(
        vcf_paths=[vcf_path], min_maf=0.0, include_missing_mafs=True
    )

    # query for all variants
    assert VALID_SIMPLE_VARIANTS_APPROX == variant_overlap_detector_query(
        variant_overlap_detector, refname="chr2", start=8000, end=9101
    )

    # query for all variants, except for the last one
    assert VALID_SIMPLE_VARIANTS_APPROX[:-1] == variant_overlap_detector_query(
        variant_overlap_detector, refname="chr2", start=8000, end=9100
    )

    # query for variants up to position 7999 (none)
    assert [] == variant_overlap_detector_query(
        variant_overlap_detector, refname="chr2", start=7000, end=7999
    )

    # query for variants between positions 8000 and 8999 (just the complex-sv variant)
    assert variant_overlap_detector_query(
        variant_overlap_detector, refname="chr2", start=8000, end=8999
    ) == get_simple_variant_approx_by_id("complex-variant-sv-1/1")

    # query for variants up to and including position 9000 (complex-sv + rare-dbsnp-snp1)
    assert variant_overlap_detector_query(
        variant_overlap_detector, refname="chr2", start=8000, end=9000
    ) == get_simple_variant_approx_by_id("complex-variant-sv-1/1", "rare-dbsnp-snp1-1/1")


@pytest.mark.parametrize("include_missing_mafs", [False, True])
def test_variant_overlap_query_maf_filter(vcf_path: Path, include_missing_mafs: bool) -> None:
    """Test that `VariantOverlapDetector.query()` MAF filtering is as expected.
    `include_missing_mafs` is parameterized in both the class constructor and in the query to
    demonstrate that it is only the query_method setting that changes the test results.
    """
    variant_overlap_detector = VariantOverlapDetector(
        vcf_paths=[vcf_path], min_maf=0.0, include_missing_mafs=include_missing_mafs
    )
    query = variant_overlap_detector_query(
        variant_overlap_detector,
        refname="chr2",
        start=8000,
        end=9101,
        maf=0.1,
        include_missing_mafs=include_missing_mafs,
    )

    if not include_missing_mafs:
        assert query == get_simple_variant_approx_by_id(
            "common-mixed-1/2",
            "common-mixed-2/2",
        )
    else:
        assert query == get_simple_variant_approx_by_id(
            "complex-variant-sv-1/1",
            "common-mixed-1/2",
            "common-mixed-2/2",
        )


@pytest.mark.parametrize("include_missing_mafs", [False, True])
def test_file_based_variant_query(vcf_path: Path, include_missing_mafs: bool) -> None:
    """Test that `FileBasedVariantLookup.query()` MAF filtering is as expected."""
    with FileBasedVariantLookup(
        vcf_paths=[vcf_path], min_maf=0.0, include_missing_mafs=include_missing_mafs
    ) as file_based_vcf_query:
        query = [
            _round_simple_variant(simple_variant)
            for simple_variant in file_based_vcf_query.query(
                refname="chr2",
                start=8000,
                end=9100,  # while "common-mixed-2/2" starts at 9101, in the VCf is starts at 9100
                maf=0.05,
                include_missing_mafs=include_missing_mafs,
            )
        ]

    if not include_missing_mafs:
        assert query == get_simple_variant_approx_by_id(
            "common-multiallelic-1/2",
            "common-multiallelic-2/2",
            "common-mixed-1/2",
            "common-mixed-2/2",
        )
    else:
        assert query == get_simple_variant_approx_by_id(
            "complex-variant-sv-1/1",
            "common-multiallelic-1/2",
            "common-multiallelic-2/2",
            "common-mixed-1/2",
            "common-mixed-2/2",
        )
