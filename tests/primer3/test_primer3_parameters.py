from dataclasses import replace

import pytest

from prymer.api.minoptmax import MinOptMax
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_parameters import Primer3Parameters


@pytest.fixture
def valid_primer3_params() -> Primer3Parameters:
    return Primer3Parameters(
        amplicon_sizes=MinOptMax(min=200, opt=250, max=300),
        amplicon_tms=MinOptMax(min=55.0, opt=60.0, max=65.0),
        primer_sizes=MinOptMax(min=18, opt=21, max=27),
        primer_tms=MinOptMax(min=55.0, opt=60.0, max=65.0),
        primer_gcs=MinOptMax(min=45.0, opt=55.0, max=60.0),
    )


def test_primer3_param_construction_valid(valid_primer3_params: Primer3Parameters) -> None:
    """Test Primer3Parameters class instantiation with valid input"""
    assert valid_primer3_params.amplicon_sizes.min == 200
    assert valid_primer3_params.amplicon_sizes.opt == 250
    assert valid_primer3_params.amplicon_sizes.max == 300
    assert valid_primer3_params.primer_gcs.min == 45.0
    assert valid_primer3_params.primer_gcs.opt == 55.0
    assert valid_primer3_params.primer_gcs.max == 60.0


def test_primer3_param_construction_raises(valid_primer3_params: Primer3Parameters) -> None:
    """Test that Primer3Parameters post_init raises with invalid input."""
    # overriding mypy here to test a case that normally would be caught by mypy
    with pytest.raises(ValueError, match="Primer Max Dinuc Bases must be an even number of bases"):
        # replace will create a new Primer instance with the provided/modified arguments
        replace(valid_primer3_params, primer_max_dinuc_bases=5)
    with pytest.raises(TypeError, match="Amplicon sizes and primer sizes must be integers"):
        replace(valid_primer3_params, amplicon_sizes=MinOptMax(min=200.0, opt=250.0, max=300.0))  # type: ignore
    with pytest.raises(TypeError, match="Amplicon sizes and primer sizes must be integers"):
        replace(valid_primer3_params, primer_sizes=MinOptMax(min=18.0, opt=21.0, max=27.0))  # type: ignore
    with pytest.raises(ValueError, match="Min primer GC-clamp must be <= max primer GC-clamp"):
        replace(valid_primer3_params, gc_clamp=(5, 0))


def test_to_input_tags_primer3_params(valid_primer3_params: Primer3Parameters) -> None:
    """Test that to_input_tags() works as expected"""
    test_dict = valid_primer3_params.to_input_tags()
    assert test_dict[Primer3InputTag.PRIMER_NUM_RETURN] == 5
    assert test_dict[Primer3InputTag.PRIMER_PRODUCT_SIZE_RANGE] == "200-300"
    assert test_dict[Primer3InputTag.PRIMER_PRODUCT_OPT_SIZE] == 250
    assert test_dict[Primer3InputTag.PRIMER_PRODUCT_MIN_TM] == 55.0
    assert test_dict[Primer3InputTag.PRIMER_PRODUCT_OPT_TM] == 60.0
    assert test_dict[Primer3InputTag.PRIMER_PRODUCT_MAX_TM] == 65.0
    assert test_dict[Primer3InputTag.PRIMER_MIN_SIZE] == 18
    assert test_dict[Primer3InputTag.PRIMER_OPT_SIZE] == 21
    assert test_dict[Primer3InputTag.PRIMER_MAX_SIZE] == 27
    assert test_dict[Primer3InputTag.PRIMER_MIN_TM] == 55.0
    assert test_dict[Primer3InputTag.PRIMER_OPT_TM] == 60.0
    assert test_dict[Primer3InputTag.PRIMER_MAX_TM] == 65.0
    assert test_dict[Primer3InputTag.PRIMER_MIN_GC] == 45.0
    assert test_dict[Primer3InputTag.PRIMER_OPT_GC_PERCENT] == 55.0
    assert test_dict[Primer3InputTag.PRIMER_MAX_GC] == 60.0
    assert test_dict[Primer3InputTag.PRIMER_GC_CLAMP] == 0
    assert test_dict[Primer3InputTag.PRIMER_MAX_END_GC] == 5
    assert test_dict[Primer3InputTag.PRIMER_MAX_POLY_X] == 5
    assert test_dict[Primer3InputTag.PRIMER_MAX_NS_ACCEPTED] == 1
    assert test_dict[Primer3InputTag.PRIMER_LOWERCASE_MASKING] == 1
    ambiguous_primer_design = replace(valid_primer3_params, avoid_masked_bases=False)
    ambiguous_dict = ambiguous_primer_design.to_input_tags()
    assert ambiguous_dict[Primer3InputTag.PRIMER_LOWERCASE_MASKING] == 0


def test_max_ampl_length(valid_primer3_params: Primer3Parameters) -> None:
    """Test that max_amplicon_length() returns expected int"""
    assert valid_primer3_params.max_amplicon_length == 300
    change_max_length = replace(
        valid_primer3_params, amplicon_sizes=MinOptMax(min=200, opt=500, max=1000)
    )
    assert change_max_length.max_amplicon_length == 1000


def test_max_primer_length(valid_primer3_params: Primer3Parameters) -> None:
    """Test that max_primer_length() returns expected int"""
    assert valid_primer3_params.max_primer_length == 27
    change_max_length = replace(
        valid_primer3_params, primer_sizes=MinOptMax(min=18, opt=35, max=50)
    )
    assert change_max_length.max_primer_length == 50
