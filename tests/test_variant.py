from pathlib import Path

import pytest
from fgpyo.vcf.builder import VariantBuilder
from fgpyo.vcf.builder import VcfFieldNumber
from fgpyo.vcf.builder import VcfFieldType
from pytest import approx

from prymer import Span
from prymer import VariantLookup
from prymer.variant import SimpleVariant


@pytest.fixture
def builder() -> VariantBuilder:
    builder = VariantBuilder(
        sd={
            "chr1": {"ID": "chr1", "length": 5000000},
            "chr2": {"ID": "chr2", "length": 5000000},
            "chr3": {"ID": "chr3", "length": 5000000},
        }
    )

    builder.add_info_header(name="CAF", field_type=VcfFieldType.STRING, number=1)
    builder.add_info_header(name="AN", field_type=VcfFieldType.INTEGER, number=1)
    builder.add_info_header(
        name="AC", field_type=VcfFieldType.INTEGER, number=VcfFieldNumber.NUM_ALT_ALLELES
    )
    builder.add_info_header(
        name="AF", field_type=VcfFieldType.FLOAT, number=VcfFieldNumber.NUM_ALT_ALLELES
    )
    builder.add_info_header(name="SVTYPE", field_type=VcfFieldType.STRING)
    builder.add_info_header(name="SVLEN", field_type=VcfFieldType.INTEGER)

    return builder


def test_simple_variant_build_simple(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, id="foo", ref="A", alts="C")
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].refname == "chr1"
    assert variants[0].start == 1000
    assert variants[0].end == 1000
    assert variants[0].ref == "A"
    assert variants[0].alt == "C"
    assert variants[0].af is None


def test_simple_variant_build_with_af(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, ref="A", alts="C", info={"AF": 0.025})
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].af == approx(0.025)

    rec = builder.add("chr1", pos=2000, ref="A", alts=["C", "G"], info={"AF": [0.1, 0.2]})
    variants = SimpleVariant.build(rec)
    assert len(variants) == 2
    assert variants[0].af == approx(0.1)
    assert variants[1].af == approx(0.2)


def test_simple_variant_build_with_caf(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, ref="A", alts="C", info={"CAF": "0.9,0.1"})
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].af == approx(0.1)

    rec = builder.add("chr1", pos=2000, ref="A", alts=["C", "G"], info={"CAF": "0.7,0.2,0.1"})
    variants = SimpleVariant.build(rec)
    assert len(variants) == 2
    assert variants[0].af == approx(0.2)
    assert variants[1].af == approx(0.1)


def test_simple_variant_build_with_ac_and_an(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, ref="A", alts="C", info={"AN": 100, "AC": [10]})
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].af == approx(0.1)

    rec = builder.add("chr1", pos=2000, ref="A", alts=["C", "G"], info={"AC": [20, 10], "AN": 100})
    variants = SimpleVariant.build(rec)
    assert len(variants) == 2
    assert variants[0].af == approx(0.2)
    assert variants[1].af == approx(0.1)


def test_simple_variant_builds_mnp(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, ref="AC", alts="CT")
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].ref == "AC"
    assert variants[0].alt == "CT"
    assert variants[0].start == 1000
    assert variants[0].end == 1001


def test_simple_variant_builds_deletion(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, ref="ACG", alts="A")
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].ref == "CG"
    assert variants[0].alt == ""
    assert variants[0].start == 1001
    assert variants[0].end == 1002


def test_simple_variant_builds_insertion(builder: VariantBuilder) -> None:
    rec = builder.add("chr1", pos=1000, ref="A", alts="ACG")
    variants = SimpleVariant.build(rec)
    assert len(variants) == 1
    assert variants[0].ref == ""
    assert variants[0].alt == "CG"
    assert variants[0].start == 1000
    assert variants[0].end == 1001


def test_variant_lookup_filtering(tmp_path: Path, builder: VariantBuilder) -> None:
    builder.add("chr1", pos=1000, ref="A", alts="C", info={"AF": 0.1})
    builder.add("chr1", pos=1010, ref="A", alts=["C", "G"], info={"AF": [0.1, 0.01]})
    builder.add("chr1", pos=1020, ref="A", alts="C")
    vcf = tmp_path / "test.vcf.gz"
    builder.to_path(vcf)

    lookup = VariantLookup([vcf], min_maf=0.1, include_missing_mafs=False)
    variants = lookup.query("chr1", 1000, 1020)
    assert len(variants) == 2
    assert variants[0].start == 1000
    assert variants[1].start == 1010 and variants[1].alt == "C"

    lookup = VariantLookup([vcf], min_maf=0.1, include_missing_mafs=True)
    variants = lookup.query("chr1", 1000, 1020)
    assert len(variants) == 3
    assert variants[0].start == 1000
    assert variants[1].start == 1010 and variants[1].alt == "C"
    assert variants[2].start == 1020

    lookup = VariantLookup([vcf], min_maf=0.001, include_missing_mafs=True)
    variants = lookup.query("chr1", 1000, 1020)
    assert len(variants) == 4


def test_variant_lookup_masking(tmp_path: Path, builder: VariantBuilder) -> None:
    builder.add("chr1", pos=1000, ref="A", alts="C", info={"AF": 0.1})
    builder.add("chr1", pos=1010, ref="ACCG", alts="A", info={"AF": 0.1})
    builder.add("chr1", pos=1020, ref="A", alts="ACGATGCA", info={"AF": 0.1})
    builder.add("chr1", pos=1030, ref="AA", alts="CG", info={"AF": 0.1})
    vcf = tmp_path / "test.vcf.gz"
    builder.to_path(vcf)

    lookup = VariantLookup([vcf], min_maf=0.1, include_missing_mafs=False)
    #              0         1         2         3         4         5
    #              012345678901234567890123456789012345678901234567890123456789
    bases_fresh = "ACGCTCGATCACCGGACTGAATGACTGACTAACTGATGCTGCGTGTCTGATGATGCTGCT"
    soft_masked = "aCGCTCGATCAccgGACTGAatGACTGACTaaCTGATGCTGCGTGTCTGATGATGCTGCT"
    hard_masked = "NCGCTCGATCANNNGACTGANNGACTGACTNNCTGATGCTGCGTGTCTGATGATGCTGCT"
    soft, hard = lookup.mask(bases=bases_fresh, region=Span("chr1", 1000, 1059))

    assert soft == soft_masked
    assert hard == hard_masked
