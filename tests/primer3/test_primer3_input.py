import pytest

from prymer.api import Span
from prymer.api import Strand
from prymer.api.minoptmax import MinOptMax
from prymer.primer3 import DesignLeftPrimersTask
from prymer.primer3 import DesignPrimerPairsTask
from prymer.primer3 import DesignRightPrimersTask
from prymer.primer3 import Primer3Input
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_parameters import PrimerAndAmpliconParameters
from prymer.primer3.primer3_parameters import ProbeParameters
from prymer.primer3.primer3_task import PickHybProbeOnly
from prymer.primer3.primer3_task import Primer3TaskType
from prymer.primer3.primer3_weights import PrimerAndAmpliconWeights
from prymer.primer3.primer3_weights import ProbeWeights


@pytest.fixture
def valid_primer_amplicon_params() -> PrimerAndAmpliconParameters:
    return PrimerAndAmpliconParameters(
        amplicon_sizes=MinOptMax(min=200, opt=250, max=300),
        amplicon_tms=MinOptMax(min=55.0, opt=60.0, max=65.0),
        primer_sizes=MinOptMax(min=18, opt=21, max=27),
        primer_tms=MinOptMax(min=55.0, opt=60.0, max=65.0),
        primer_gcs=MinOptMax(min=45.0, opt=55.0, max=60.0),
    )


@pytest.fixture
def valid_probe_params() -> ProbeParameters:
    return ProbeParameters(
        probe_sizes=MinOptMax(min=18, opt=22, max=30),
        probe_tms=MinOptMax(min=65.0, opt=70.0, max=75.0),
        probe_gcs=MinOptMax(min=45.0, opt=55.0, max=60.0),
    )


@pytest.fixture
def valid_primer_weights() -> PrimerAndAmpliconWeights:
    return PrimerAndAmpliconWeights()


@pytest.fixture
def valid_probe_weights() -> ProbeWeights:
    return ProbeWeights(probe_wt_hairpin_th=100.0)


@pytest.mark.parametrize(
    "task_type",
    [
        DesignRightPrimersTask(),
        DesignLeftPrimersTask(),
        DesignPrimerPairsTask(),
    ],
)
def test_primer_design_only_valid(
    valid_primer_amplicon_params: PrimerAndAmpliconParameters,
    task_type: Primer3TaskType,
    valid_primer_weights: PrimerAndAmpliconWeights,
) -> None:
    test_design_region = Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE)
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    test_input = Primer3Input(
        target=test_target,
        primer_weights=valid_primer_weights,
        task=task_type,
        primer_and_amplicon_params=valid_primer_amplicon_params,
    )
    mapped_dict = test_input.to_input_tags(design_region=test_design_region)
    assert len(mapped_dict.keys()) == 46


@pytest.mark.parametrize(
    "task_type", [DesignRightPrimersTask(), DesignLeftPrimersTask(), DesignPrimerPairsTask()]
)
def test_primer_design_only_raises(
    task_type: Primer3TaskType, valid_primer_weights: PrimerAndAmpliconWeights
) -> None:
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    with pytest.raises(ValueError, match="Primer3 requires at least one set of parameters"):
        Primer3Input(
            target=test_target,
            primer_weights=valid_primer_weights,
            task=task_type,
            primer_and_amplicon_params=None,
        )


def test_probe_design_only_valid(
    valid_probe_params: ProbeParameters, valid_probe_weights: ProbeWeights
) -> None:
    test_design_region = Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE)
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    test_input = Primer3Input(
        target=test_target,
        probe_weights=valid_probe_weights,
        task=PickHybProbeOnly(),
        probe_params=valid_probe_params,
        primer_and_amplicon_params=None,
    )
    mapped_dict = test_input.to_input_tags(design_region=test_design_region)
    assert mapped_dict[Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO] == 1
    assert mapped_dict[Primer3InputTag.PRIMER_INTERNAL_WT_HAIRPIN_TH] == 100.0

    assert len(mapped_dict.keys()) == 29

    # test instantiation of default `ProbeWeights` when they are not provided
    altered_input = Primer3Input(
        target=test_target,
        probe_weights=None,
        task=PickHybProbeOnly(),
        probe_params=valid_probe_params,
        primer_and_amplicon_params=None,
    )
    altered_mapped_dict = altered_input.to_input_tags(design_region=test_target)
    assert altered_mapped_dict[Primer3InputTag.PRIMER_INTERNAL_WT_HAIRPIN_TH] == 1.0


def test_probe_design_only_raises(valid_probe_weights: ProbeWeights) -> None:
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    with pytest.raises(ValueError, match="Primer3 requires at least one set"):
        Primer3Input(
            target=test_target,
            probe_weights=valid_probe_weights,
            task=PickHybProbeOnly(),
            primer_and_amplicon_params=None,
        )


@pytest.mark.parametrize(
    "task_type",
    [
        DesignRightPrimersTask(),
        DesignLeftPrimersTask(),
        DesignPrimerPairsTask(),
        PickHybProbeOnly(),
    ],
)
def test_no_params_given_raises(
    valid_primer_weights: PrimerAndAmpliconWeights, task_type: Primer3TaskType
) -> None:
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    with pytest.raises(ValueError, match="Primer3 requires at least one set"):
        Primer3Input(
            target=test_target,
            primer_weights=valid_primer_weights,
            task=task_type,
            primer_and_amplicon_params=None,
            probe_params=None,
        )


@pytest.mark.parametrize(
    "task_type, expected_req_primer_amplicon_params, expected_req_probe_params",
    [
        (DesignPrimerPairsTask(), True, False),
        (DesignRightPrimersTask(), True, False),
        (DesignLeftPrimersTask(), True, False),
        (PickHybProbeOnly(), False, True),
    ],
)
def test_requires_params_sets(
    task_type: Primer3TaskType,
    valid_probe_params: ProbeParameters,
    valid_primer_amplicon_params: PrimerAndAmpliconParameters,
    valid_primer_weights: PrimerAndAmpliconWeights,
    valid_probe_weights: ProbeWeights,
    expected_req_primer_amplicon_params: bool,
    expected_req_probe_params: bool,
) -> None:
    test_target = Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE)
    test_input = Primer3Input(
        target=test_target,
        primer_weights=valid_primer_weights,
        probe_weights=valid_probe_weights,
        task=task_type,
        probe_params=valid_probe_params,
        primer_and_amplicon_params=valid_primer_amplicon_params,
    )
    assert test_input.task.requires_probe_params == expected_req_probe_params
    assert test_input.task.requires_primer_amplicon_params == expected_req_primer_amplicon_params
