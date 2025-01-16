from dataclasses import replace

import pytest

from prymer import MinOptMax
from prymer import Span
from prymer import Strand
from prymer.primer3 import DesignLeftPrimersTask
from prymer.primer3 import DesignPrimerPairsTask
from prymer.primer3 import DesignRightPrimersTask
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_parameters import AmpliconParameters
from prymer.primer3.primer3_parameters import ProbeParameters
from prymer.primer3.primer3_task import PickHybProbeOnly
from prymer.primer3.primer3_task import Primer3TaskType


@pytest.fixture
def target() -> Span:
    return Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)


@pytest.fixture
def design_region() -> Span:
    return Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE)


@pytest.fixture
def valid_primer_amplicon_params(target: Span) -> AmpliconParameters:
    return AmpliconParameters(
        target=target,
        task=DesignLeftPrimersTask(),
        amplicon_sizes=MinOptMax(min=200, opt=250, max=300),
        amplicon_tms=MinOptMax(min=55.0, opt=60.0, max=65.0),
        primer_sizes=MinOptMax(min=18, opt=21, max=27),
        primer_tms=MinOptMax(min=55.0, opt=60.0, max=65.0),
        primer_gcs=MinOptMax(min=45.0, opt=55.0, max=60.0),
    )


@pytest.fixture
def valid_probe_params(target: Span) -> ProbeParameters:
    return ProbeParameters(
        target=target,
        task=PickHybProbeOnly(),
        probe_sizes=MinOptMax(min=18, opt=22, max=30),
        probe_tms=MinOptMax(min=65.0, opt=70.0, max=75.0),
        probe_gcs=MinOptMax(min=45.0, opt=55.0, max=60.0),
    )


@pytest.mark.parametrize(
    "task_type",
    [
        DesignRightPrimersTask(),
        DesignLeftPrimersTask(),
        DesignPrimerPairsTask(),
    ],
)
def test_primer_design_only_valid(
    valid_primer_amplicon_params: AmpliconParameters,
    task_type: Primer3TaskType,
) -> None:
    test_design_region = Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE)
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    test_input = replace(
        valid_primer_amplicon_params,
        target=test_target,
        task=task_type,
    )
    mapped_dict = test_input.to_input_tags(design_region=test_design_region)
    assert len(mapped_dict.keys()) == 44


def test_probe_design_only_valid(
    valid_probe_params: ProbeParameters,
) -> None:
    test_design_region = Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE)
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    test_input = replace(
        valid_probe_params,
        target=test_target,
        task=PickHybProbeOnly(),
    )
    mapped_dict = test_input.to_input_tags(design_region=test_design_region)
    assert mapped_dict[Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO] == 1
    assert Primer3InputTag.PRIMER_NUM_RETURN in mapped_dict

    assert len(mapped_dict.keys()) == 28

    # test instantiation of default `ProbeParameters` when they are not provided
    altered_input = replace(
        valid_probe_params,
        target=test_target,
        task=PickHybProbeOnly(),
    )
    altered_mapped_dict = altered_input.to_input_tags(design_region=test_target)
    assert altered_mapped_dict[Primer3InputTag.PRIMER_INTERNAL_WT_GC_PERCENT_GT] == 0.0


def test_primer_amplicon_param_construction_valid(
    valid_primer_amplicon_params: AmpliconParameters,
) -> None:
    """Test PrimerAndAmpliconParameters class instantiation with valid input"""
    assert valid_primer_amplicon_params.amplicon_sizes.min == 200
    assert valid_primer_amplicon_params.amplicon_sizes.opt == 250
    assert valid_primer_amplicon_params.amplicon_sizes.max == 300
    assert valid_primer_amplicon_params.primer_gcs.min == 45.0
    assert valid_primer_amplicon_params.primer_gcs.opt == 55.0
    assert valid_primer_amplicon_params.primer_gcs.max == 60.0


def test_probe_param_construction_valid(
    valid_probe_params: ProbeParameters,
) -> None:
    """Test ProbeParameters class instantiation with valid input. Assert that the *_thermo fields
    (which were not explicitly set) are set to `valid_probe_params.probe_tms.min` - 10.0."""
    expected_thermo_value: float = valid_probe_params.probe_tms.min - 10.0
    assert valid_probe_params.probe_sizes.min == 18
    assert valid_probe_params.probe_sizes.opt == 22
    assert valid_probe_params.probe_sizes.max == 30
    assert valid_probe_params.probe_tms.min == 65.0
    assert valid_probe_params.probe_tms.opt == 70.0
    assert valid_probe_params.probe_tms.max == 75.0
    assert valid_probe_params.probe_gcs.min == 45.0
    assert valid_probe_params.probe_gcs.opt == 55.0
    assert valid_probe_params.probe_gcs.max == 60.0
    assert valid_probe_params.probe_max_homodimer_tm == expected_thermo_value
    assert valid_probe_params.probe_max_3p_homodimer_tm == expected_thermo_value
    assert valid_probe_params.probe_max_hairpin_tm == expected_thermo_value


def test_primer_amplicon_param_construction_raises(
    valid_primer_amplicon_params: AmpliconParameters,
) -> None:
    """Test that PrimerAndAmpliconParameters post_init raises with invalid input."""
    # overriding mypy here to test a case that normally would be caught by mypy
    with pytest.raises(ValueError, match="Primer Max Dinuc Bases must be an even number of bases"):
        # replace will create a new Primer instance with the provided/modified arguments
        replace(valid_primer_amplicon_params, primer_max_dinuc_bases=5)
    with pytest.raises(TypeError, match="Amplicon sizes and primer sizes must be integers"):
        replace(
            valid_primer_amplicon_params,
            amplicon_sizes=MinOptMax(min=200.0, opt=250.0, max=300.0),  # type: ignore
        )
    with pytest.raises(TypeError, match="Amplicon sizes and primer sizes must be integers"):
        replace(valid_primer_amplicon_params, primer_sizes=MinOptMax(min=18.0, opt=21.0, max=27.0))  # type: ignore
    with pytest.raises(ValueError, match="Min primer GC-clamp must be <= max primer GC-clamp"):
        replace(valid_primer_amplicon_params, gc_clamp=(5, 0))


def test_primer_probe_param_construction_raises(
    valid_probe_params: ProbeParameters,
) -> None:
    """Test that Primer3Parameters post_init raises with invalid input."""
    # overriding mypy here to test a case that normally would be caught by mypy
    with pytest.raises(TypeError, match="Probe sizes must be integers"):
        # replace will create a new Primer instance with the provided/modified arguments
        # we use `type: ignore` here to bypass mypy
        replace(
            valid_probe_params,
            probe_sizes=MinOptMax(min=18.1, opt=22.1, max=30.1),  # type: ignore
        )
    with pytest.raises(TypeError, match="Probe melting temperatures and GC content must be floats"):
        replace(valid_probe_params, probe_tms=MinOptMax(min=55, opt=60, max=65))


def test_primer_amplicon_params_to_input_tags(
    valid_primer_amplicon_params: AmpliconParameters,
    design_region: Span,
) -> None:
    """Test that to_input_tags() works as expected"""
    test_dict = valid_primer_amplicon_params.to_input_tags(design_region=design_region)
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
    ambiguous_primer_design = replace(valid_primer_amplicon_params, avoid_masked_bases=False)
    ambiguous_dict = ambiguous_primer_design.to_input_tags(design_region=design_region)
    assert ambiguous_dict[Primer3InputTag.PRIMER_LOWERCASE_MASKING] == 0


def test_max_ampl_length(valid_primer_amplicon_params: AmpliconParameters) -> None:
    """Test that max_amplicon_length() returns expected int"""
    assert valid_primer_amplicon_params.max_amplicon_length == 300
    change_max_length = replace(
        valid_primer_amplicon_params, amplicon_sizes=MinOptMax(min=200, opt=500, max=1000)
    )
    assert change_max_length.max_amplicon_length == 1000


def test_max_primer_length(valid_primer_amplicon_params: AmpliconParameters) -> None:
    """Test that max_primer_length() returns expected int"""
    assert valid_primer_amplicon_params.max_primer_length == 27
    change_max_length = replace(
        valid_primer_amplicon_params, primer_sizes=MinOptMax(min=18, opt=35, max=50)
    )
    assert change_max_length.max_primer_length == 50


def test_primer_weights_valid(
    design_region: Span, valid_primer_amplicon_params: AmpliconParameters
) -> None:
    """Test instantiation of `AmpliconParameters` object with valid input"""
    test_dict = valid_primer_amplicon_params.to_input_tags(design_region=design_region)
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT] == 1
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_GT] == 1
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_TM_LT] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_TM_GT] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_WT_END_STABILITY] == 0.25
    assert test_dict[Primer3InputTag.PRIMER_WT_GC_PERCENT_LT] == 0.25
    assert test_dict[Primer3InputTag.PRIMER_WT_GC_PERCENT_GT] == 0.25
    assert test_dict[Primer3InputTag.PRIMER_WT_SELF_ANY] == 0.1
    assert test_dict[Primer3InputTag.PRIMER_WT_SELF_END] == 0.1
    assert test_dict[Primer3InputTag.PRIMER_WT_SIZE_LT] == 0.5
    assert test_dict[Primer3InputTag.PRIMER_WT_SIZE_GT] == 0.1
    assert test_dict[Primer3InputTag.PRIMER_WT_TM_LT] == 1.0
    assert test_dict[Primer3InputTag.PRIMER_WT_TM_GT] == 1.0
    assert len((test_dict.values())) == 44


def test_probe_weights_valid(design_region: Span, valid_probe_params: ProbeParameters) -> None:
    test_dict = valid_probe_params.to_input_tags(design_region=design_region)
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_SIZE_LT] == 1.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_SIZE_GT] == 1.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_TM_LT] == 1.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_TM_GT] == 1.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_GC_PERCENT_LT] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_GC_PERCENT_GT] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_SELF_ANY_TH] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_SELF_END_TH] == 0.0
    assert test_dict[Primer3InputTag.PRIMER_INTERNAL_WT_HAIRPIN_TH] == 0.0
    assert len(test_dict) == 28


def test_primer_weights_to_input_tags(
    design_region: Span, valid_primer_amplicon_params: AmpliconParameters
) -> None:
    """Test results from to_input_tags() with and without default values"""
    default_map = valid_primer_amplicon_params.to_input_tags(design_region)
    assert default_map[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT] == 1
    customized_map = replace(valid_primer_amplicon_params, product_size_lt=5).to_input_tags(
        design_region=design_region
    )
    assert customized_map[Primer3InputTag.PRIMER_PAIR_WT_PRODUCT_SIZE_LT] == 5
