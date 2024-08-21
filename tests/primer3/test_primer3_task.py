import pytest

from prymer.api.span import Span
from prymer.api.span import Strand
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_task import DesignLeftPrimersTask
from prymer.primer3.primer3_task import DesignPrimerPairsTask
from prymer.primer3.primer3_task import DesignRightPrimersTask
from prymer.primer3.primer3_task import Primer3Task


def test_pair_design_construction() -> None:
    test_pair_design = DesignPrimerPairsTask()
    assert test_pair_design.task_type == "PAIR"
    assert test_pair_design.count_tag == "PRIMER_PAIR_NUM_RETURNED"


def test_right_design_construction() -> None:
    test_right_primer_design = DesignRightPrimersTask()
    assert test_right_primer_design.task_type == "RIGHT"
    assert test_right_primer_design.count_tag == "PRIMER_RIGHT_NUM_RETURNED"


def test_left_design_construction() -> None:
    test_left_primer_design = DesignLeftPrimersTask()
    assert test_left_primer_design.task_type == "LEFT"
    assert test_left_primer_design.count_tag == "PRIMER_LEFT_NUM_RETURNED"


@pytest.mark.parametrize(
    "design_region, target_region, expected_seq_target",
    [
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE),
            "200,101",
        ),
        (
            Span(refname="chr1", start=200, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=220, end=480, strand=Strand.POSITIVE),
            "21,261",
        ),
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=1, end=500, strand=Strand.NEGATIVE),
            "1,500",
        ),
    ],
)
def test_pair_design_to_input_tags(
    design_region: Span, target_region: Span, expected_seq_target: str
) -> None:
    test_task = DesignPrimerPairsTask()
    test_tags = test_task.to_input_tags(design_region=design_region, target=target_region)
    assert test_tags[Primer3InputTag.SEQUENCE_TARGET] == expected_seq_target
    assert test_tags[Primer3InputTag.PRIMER_TASK] == "generic"
    assert test_tags[Primer3InputTag.PRIMER_PICK_LEFT_PRIMER] == 1
    assert test_tags[Primer3InputTag.PRIMER_PICK_RIGHT_PRIMER] == 1
    assert test_tags[Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO] == 0


@pytest.mark.parametrize(
    "design_region, target_region, expected_seq_target",
    [
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE),
            "1,199",
        ),
        (
            Span(refname="chr1", start=200, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=220, end=480, strand=Strand.POSITIVE),
            "1,20",
        ),
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=1, end=500, strand=Strand.NEGATIVE),
            "1,0",
        ),
    ],
)
def test_left_design_to_input_tags(
    design_region: Span, target_region: Span, expected_seq_target: str
) -> None:
    test_task = DesignLeftPrimersTask()
    test_tags = test_task.to_input_tags(design_region=design_region, target=target_region)
    assert test_tags[Primer3InputTag.SEQUENCE_INCLUDED_REGION] == expected_seq_target
    assert test_tags[Primer3InputTag.PRIMER_TASK] == "pick_primer_list"
    assert test_tags[Primer3InputTag.PRIMER_PICK_LEFT_PRIMER] == 1
    assert test_tags[Primer3InputTag.PRIMER_PICK_RIGHT_PRIMER] == 0
    assert test_tags[Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO] == 0


@pytest.mark.parametrize(
    "design_region, target_region, expected_seq_target",
    [
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE),
            "300,200",
        ),
        (
            Span(refname="chr1", start=200, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=220, end=480, strand=Strand.POSITIVE),
            "281,20",
        ),
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=1, end=500, strand=Strand.NEGATIVE),
            "500,0",
        ),
    ],
)
def test_right_design_to_input_tags(
    design_region: Span, target_region: Span, expected_seq_target: str
) -> None:
    test_task = DesignRightPrimersTask()
    test_tags = test_task.to_input_tags(design_region=design_region, target=target_region)
    assert test_tags[Primer3InputTag.SEQUENCE_INCLUDED_REGION] == expected_seq_target
    assert test_tags[Primer3InputTag.PRIMER_TASK] == "pick_primer_list"
    assert test_tags[Primer3InputTag.PRIMER_PICK_LEFT_PRIMER] == 0
    assert test_tags[Primer3InputTag.PRIMER_PICK_RIGHT_PRIMER] == 1
    assert test_tags[Primer3InputTag.PRIMER_PICK_INTERNAL_OLIGO] == 0


@pytest.mark.parametrize(
    "design_region, target_region",
    [
        (
            Span(refname="chr1", start=1, end=500, strand=Strand.POSITIVE),
            Span(refname="chr1", start=400, end=600, strand=Strand.POSITIVE),
        ),  # target region is outside design region, right side
        (
            Span(refname="chr1", start=1, end=10, strand=Strand.POSITIVE),
            Span(refname="chr1", start=9, end=20, strand=Strand.POSITIVE),
        ),  # target region outside design region, right side
        (
            Span(refname="chr1", start=100, end=200, strand=Strand.POSITIVE),
            Span(refname="chr1", start=200, end=300, strand=Strand.POSITIVE),
        ),  # target region outside design region, right side
        (
            Span(refname="chr1", start=100, end=200, strand=Strand.POSITIVE),
            Span(refname="chr1", start=1, end=101, strand=Strand.POSITIVE),
        ),  # target region outside design, left side
        (
            Span(refname="chr1", start=100, end=200, strand=Strand.POSITIVE),
            Span(refname="chr1", start=1, end=100, strand=Strand.POSITIVE),
        ),  # target region outside design, left side
    ],
)
@pytest.mark.parametrize(
    "task_type", [DesignRightPrimersTask(), DesignLeftPrimersTask(), DesignPrimerPairsTask()]
)
def test_invalid_design_target_combo(
    design_region: Span,
    target_region: Span,
    task_type: Primer3Task,
) -> None:
    """Test that all Primer3Tasks raise an error when the target region is not wholly
    contained by the design region."""
    assert design_region.overlaps(target_region)
    test_task = task_type
    with pytest.raises(ValueError, match="Target not contained within design region:"):
        test_task.to_input_tags(design_region=design_region, target=target_region)
