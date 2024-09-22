import logging
from dataclasses import replace
from pathlib import Path

import pysam
import pytest
from fgpyo.sequence import reverse_complement

from prymer.api.minoptmax import MinOptMax
from prymer.api.primer import Primer
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.api.span import Strand
from prymer.api.variant_lookup import cached
from prymer.primer3.primer3 import Primer3
from prymer.primer3.primer3 import Primer3Failure
from prymer.primer3.primer3 import Primer3Result
from prymer.primer3.primer3 import _has_acceptable_dinuc_run
from prymer.primer3.primer3_input import Primer3Input
from prymer.primer3.primer3_parameters import PrimerAndAmpliconParameters
from prymer.primer3.primer3_parameters import ProbeParameters
from prymer.primer3.primer3_task import DesignLeftPrimersTask
from prymer.primer3.primer3_task import DesignPrimerPairsTask
from prymer.primer3.primer3_task import DesignRightPrimersTask
from prymer.primer3.primer3_task import PickHybProbeOnly


@pytest.fixture(scope="session")
def genome_ref() -> Path:
    return Path(__file__).parent / "data" / "miniref.fa"


@pytest.fixture
def vcf_path() -> Path:
    return Path(__file__).parent / "data" / "miniref.variants.vcf.gz"


@pytest.fixture
def valid_probe_params_no_exclude() -> ProbeParameters:
    return ProbeParameters(
        probe_sizes=MinOptMax(min=18, max=30, opt=22),
        probe_tms=MinOptMax(min=55.0, max=100.0, opt=70.0),
        probe_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0),
    )


@pytest.fixture
def single_primer_params() -> PrimerAndAmpliconParameters:
    return PrimerAndAmpliconParameters(
        amplicon_sizes=MinOptMax(min=100, max=250, opt=200),
        amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=70.0),
        primer_sizes=MinOptMax(min=29, max=31, opt=30),
        primer_tms=MinOptMax(min=63.0, max=67.0, opt=65.0),
        primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0),
        primer_max_polyX=4,
        number_primers_return=1000,
    )


@pytest.fixture
def pair_primer_params() -> PrimerAndAmpliconParameters:
    return PrimerAndAmpliconParameters(
        amplicon_sizes=MinOptMax(min=100, max=200, opt=150),
        amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=72.5),
        primer_sizes=MinOptMax(min=20, max=30, opt=25),
        primer_tms=MinOptMax(min=55.0, max=75.0, opt=65.0),
        primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0),
        primer_max_polyX=4,
        number_primers_return=10,
    )


@pytest.fixture
def design_fail_gen_primer3_params() -> PrimerAndAmpliconParameters:
    return PrimerAndAmpliconParameters(
        amplicon_sizes=MinOptMax(min=200, max=300, opt=250),
        amplicon_tms=MinOptMax(min=65.0, max=75.0, opt=74.0),
        primer_sizes=MinOptMax(min=24, max=27, opt=26),
        primer_tms=MinOptMax(min=65.0, max=75.0, opt=74.0),
        primer_gcs=MinOptMax(min=55.0, max=65.0, opt=62.0),
    )


def make_primer(bases: str, refname: str, start: int, end: int) -> Primer:
    return Primer(
        bases=bases,
        tm=55,
        penalty=5,
        span=Span(refname=refname, start=start, end=end),
    )


def make_primer_pair(left: Primer, right: Primer, genome_ref: Path) -> PrimerPair:
    ref = pysam.FastaFile(str(genome_ref))  # pysam expects a str instead of Path
    amplicon_span = Span(
        refname=left.span.refname,
        start=left.span.start,
        end=right.span.end,
    )
    amplicon_sequence = ref.fetch(
        region=f"{amplicon_span.refname}:{amplicon_span.start}-{amplicon_span.end}"
    )
    return PrimerPair(
        left_primer=left,
        right_primer=right,
        amplicon_sequence=amplicon_sequence,
        penalty=5.0,
        amplicon_tm=60.0,
    )


@pytest.fixture(scope="session")
def valid_left_primers() -> list[Primer]:
    lefts: list[Primer] = [
        make_primer(bases="ACATTTGCTTCTGACACAAC", refname="chr1", start=1, end=20),
        make_primer(bases="TGTGTTCACTAGCAACCTCA", refname="chr1", start=21, end=40),
    ]
    return lefts


@pytest.fixture(scope="session")
def valid_right_primers() -> list[Primer]:
    rights: list[Primer] = [
        make_primer(
            bases=reverse_complement("TCAAGGTTACAAGACAGGTT"), refname="chr1", start=150, end=169
        ),
        make_primer(
            bases=reverse_complement("TAAGGAGACCAATAGAAACT"), refname="chr1", start=170, end=189
        ),
    ]
    return rights


@pytest.fixture(scope="session")
def valid_primer_pairs(
    valid_left_primers: list[Primer], valid_right_primers: list[Primer], genome_ref: Path
) -> list[PrimerPair]:
    primer_pairs = [
        make_primer_pair(left=left, right=right, genome_ref=genome_ref)
        for left, right in list(zip(valid_left_primers, valid_right_primers, strict=True))
    ]
    return primer_pairs


def test_design_primers_raises(
    genome_ref: Path,
    single_primer_params: PrimerAndAmpliconParameters,
) -> None:
    """Test that design_primers() raises when given an invalid argument."""

    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)

    illegal_primer3_params = replace(
        single_primer_params,
        number_primers_return="invalid",  # type: ignore
    )
    invalid_design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=illegal_primer3_params,
        task=DesignLeftPrimersTask(),
    )
    with pytest.raises(ValueError, match="Primer3 failed"):
        Primer3(genome_fasta=genome_ref).design_oligos(design_input=invalid_design_input)
    # TODO: add other Value Errors


def test_internal_probe_valid_designs(
    genome_ref: Path,
    valid_probe_params_no_exclude: ProbeParameters,
) -> None:
    """Test that internal probe designs are within the specified design specifications."""
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
    assert valid_probe_params_no_exclude is not None
    design_input = Primer3Input(
        target=target,
        probe_params=valid_probe_params_no_exclude,
        task=PickHybProbeOnly(),
    )
    with Primer3(genome_fasta=genome_ref) as designer:
        primer3_result = designer.design_oligos(design_input=design_input)
    assert len(primer3_result.filtered_designs) == 5
    for probe_design in primer3_result.filtered_designs:
        assert probe_design.self_any_th < valid_probe_params_no_exclude.probe_max_self_any_thermo
        assert probe_design.self_end_th < valid_probe_params_no_exclude.probe_max_self_end_thermo
        assert probe_design.hairpin_th < valid_probe_params_no_exclude.probe_max_hairpin_thermo
        assert (
            probe_design.longest_dinucleotide_run_length()
            <= valid_probe_params_no_exclude.probe_max_dinuc_bases
        )
        assert probe_design.span.start >= target.start
        assert probe_design.span.end <= target.end


def test_left_primer_valid_designs(
    genome_ref: Path,
    single_primer_params: PrimerAndAmpliconParameters,
) -> None:
    """Test that left primer designs are within the specified design specifications."""
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
    design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=single_primer_params,
        task=DesignLeftPrimersTask(),
    )

    with Primer3(genome_fasta=genome_ref) as designer:
        for _ in range(10):  # run many times to ensure we can re-use primer3
            left_result = designer.design_oligos(design_input=design_input)
            designed_lefts: list[Primer] = left_result.primers()
            assert all(isinstance(design, Primer) for design in designed_lefts)
            for actual_design in designed_lefts:
                assert (
                    actual_design.longest_dinucleotide_run_length()
                    <= single_primer_params.primer_max_dinuc_bases
                )
                assert (
                    single_primer_params.primer_sizes.min
                    <= actual_design.length
                    <= single_primer_params.primer_sizes.max
                )
                assert (
                    single_primer_params.primer_tms.min
                    <= actual_design.tm
                    <= single_primer_params.primer_tms.max
                )
                assert (
                    single_primer_params.primer_gcs.min
                    <= actual_design.percent_gc_content
                    <= single_primer_params.primer_gcs.max
                )
                assert actual_design.span.start < actual_design.span.end
                assert actual_design.span.end < target.start
                underlying_ref_seq = designer._fasta.fetch(  # pysam is 0-based, half-open
                    reference=actual_design.span.refname,
                    start=actual_design.span.start - 1,
                    end=actual_design.span.end,
                )
                assert actual_design.bases == underlying_ref_seq
        assert designer.is_alive


def test_right_primer_valid_designs(
    genome_ref: Path,
    single_primer_params: PrimerAndAmpliconParameters,
) -> None:
    """Test that right primer designs are within the specified design specifications."""
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
    design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=single_primer_params,
        task=DesignRightPrimersTask(),
    )
    with Primer3(genome_fasta=genome_ref) as designer:
        for _ in range(10):  # run many times to ensure we can re-use primer3
            right_result: Primer3Result = designer.design_oligos(design_input=design_input)
            designed_rights: list[Primer] = right_result.primers()
            assert all(isinstance(design, Primer) for design in designed_rights)

            for actual_design in designed_rights:
                assert (
                    actual_design.longest_dinucleotide_run_length()
                    <= single_primer_params.primer_max_dinuc_bases
                )
                assert (
                    single_primer_params.primer_sizes.min
                    <= actual_design.length
                    <= single_primer_params.primer_sizes.max
                )
                assert (
                    single_primer_params.primer_tms.min
                    <= actual_design.tm
                    <= single_primer_params.primer_tms.max
                )
                assert (
                    single_primer_params.primer_gcs.min
                    <= actual_design.percent_gc_content
                    <= single_primer_params.primer_gcs.max
                )
                assert actual_design.span.start < actual_design.span.end
                assert actual_design.span.end > actual_design.span.start
                assert actual_design.span.start > 250
                underlying_ref_seq = designer._fasta.fetch(  # pysam is 0-based, half-open
                    reference=actual_design.span.refname,
                    start=actual_design.span.start - 1,
                    end=actual_design.span.end,
                )
                assert actual_design.bases == reverse_complement(underlying_ref_seq)
        assert designer.is_alive


def test_primer_pair_design(
    genome_ref: Path, pair_primer_params: PrimerAndAmpliconParameters
) -> None:
    """Test that paired primer design produces left and right primers within design constraints.
    Additionally, assert that `PrimerPair.amplicon_sequence()` matches reference sequence."""
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
    design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=pair_primer_params,
        task=DesignPrimerPairsTask(),
    )
    with Primer3(genome_fasta=genome_ref) as designer:
        pair_result: Primer3Result = designer.design_oligos(design_input=design_input)
        designed_pairs: list[PrimerPair] = pair_result.primer_pairs()
        assert all(isinstance(design, PrimerPair) for design in designed_pairs)
        lefts = [primer_pair.left_primer for primer_pair in designed_pairs]
        rights = [primer_pair.right_primer for primer_pair in designed_pairs]

        assert len(lefts) == 10
        assert len(rights) == 10
        for pair_design in designed_pairs:
            if pair_design.amplicon_sequence is not None:
                assert len(pair_design.amplicon_sequence) <= pair_primer_params.amplicon_sizes.max
                assert (
                    pair_design.amplicon_sequence.upper()
                    == designer._fasta.fetch(  # pysam is 0-based, half-open
                        reference=pair_design.left_primer.span.refname,
                        start=pair_design.left_primer.span.start - 1,
                        end=pair_design.right_primer.span.end,
                    ).upper()
                )
            # check left primers
            assert (
                pair_primer_params.primer_sizes.min
                <= pair_design.left_primer.length
                <= pair_primer_params.primer_sizes.max
            )
            assert (
                pair_primer_params.primer_tms.min
                <= pair_design.left_primer.tm
                <= pair_primer_params.primer_tms.max
            )
            assert (
                pair_primer_params.primer_gcs.min
                <= pair_design.left_primer.percent_gc_content
                <= pair_primer_params.primer_gcs.max
            )
            assert pair_design.left_primer.bases is not None
            # check right primers
            assert (
                pair_primer_params.primer_sizes.min
                <= pair_design.right_primer.length
                <= pair_primer_params.primer_sizes.max
            )
            assert (
                pair_primer_params.primer_tms.min
                <= pair_design.right_primer.tm
                <= pair_primer_params.primer_tms.max
            )
            assert (
                pair_primer_params.primer_gcs.min
                <= pair_design.right_primer.percent_gc_content
                <= pair_primer_params.primer_gcs.max
            )
            assert pair_design.right_primer.bases is not None

            left_from_ref = designer._fasta.fetch(  # pysam is 0-based, half-open
                reference=pair_design.left_primer.span.refname,
                start=pair_design.left_primer.span.start - 1,
                end=pair_design.left_primer.span.end,
            )

            right_from_ref = reverse_complement(
                designer._fasta.fetch(  # pysam is 0-based, half-open
                    reference=pair_design.right_primer.span.refname,
                    start=pair_design.right_primer.span.start - 1,
                    end=pair_design.right_primer.span.end,
                )
            )
            assert pair_design.left_primer.bases.upper() == left_from_ref.upper()
            assert pair_design.right_primer.bases.upper() == right_from_ref.upper()


def test_fasta_close_valid(
    genome_ref: Path, single_primer_params: PrimerAndAmpliconParameters
) -> None:
    """Test that fasta file is closed when underlying subprocess is terminated."""
    designer = Primer3(genome_fasta=genome_ref)
    assert designer._fasta.is_open()
    designer.close()
    assert designer._fasta.closed
    assert not designer.is_alive
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
    design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=single_primer_params,
        task=DesignLeftPrimersTask(),
    )

    with pytest.raises(
        RuntimeError, match="Error, trying to use a subprocess that has already been terminated"
    ):
        designer.design_oligos(design_input=design_input)


@pytest.mark.parametrize(
    "region, expected_hard_masked, expected_soft_masked",
    [
        (
            Span(refname="chr2", start=9000, end=9110),
            # 9000      9010      9020      9030      9040      9050      9060      9070      9080      9090      9100      9110 # noqa
            "AATATTCTTGNTGCTTATGCNGCTGACATTGTTGCCCTCCCTAAAGCAACNAAGTAGCCTNTATTTCCCANAGTGAAAGANNACGCTGGCNNNTCAGTTANNNTACAAAAG",
            "AATATTCTTGCTGCTTATGCAGCTGACATTGTTGCCCTCCCTAAAGCAACCAAGTAGCCTTTATTTCCCACAGTGAAAGAAAACGCTGGCCTATCAGTTACATTACAAAAG",
        ),  # expected masked positions: 9010, 9020, 9050, 9060, 9070,
        # 9080 (2bp insertion: 3 bases), 9090 (2bp deletion: 2 bases), 9100 (mixed: 3 bases)
        # do not expect positions 9000 (MAF = 0.001), 9030 (MAF = 0.001), or 9040 (MAF = 0.0004814)
        # to be masked  (MAF below the provided min_maf)
        (
            Span(refname="chr2", start=9095, end=9120),
            "AGTTANNNTACAAAAGGCAGATTTCA",
            "AGTTACATTACAAAAGGCAGATTTCA",
        ),
        # 9100 (common-mixed -- alt1: CA->GG, and alt2: CA->CACACA).  The first alt masks the
        # positions [9100,9101], and the second alt masks the positions [9100,9102] (an extra
        # base for the insertion).  But the second alt is not added to variant lookup, while the
        # first variant is classified as OTHER, so [9100,9102] are masked.  FIXME: this could be
        # improved by more faithfully parsing the input VCF and representing each alternate as its
        # own simple variant.
    ],
)
def test_variant_lookup(
    genome_ref: Path,
    vcf_path: Path,
    region: Span,
    expected_hard_masked: str,
    expected_soft_masked: str,
) -> None:
    """Test that MAF filtering and masking are working as expected."""
    with Primer3(
        genome_fasta=genome_ref, variant_lookup=cached([vcf_path], min_maf=0.01)
    ) as designer:
        actual_soft_masked, actual_hard_masked = designer.get_design_sequences(region=region)
    assert actual_hard_masked == expected_hard_masked
    assert actual_soft_masked == expected_soft_masked

    # with no variant lookup should all be soft-masked
    with Primer3(genome_fasta=genome_ref, variant_lookup=None) as designer:
        actual_soft_masked, actual_hard_masked = designer.get_design_sequences(region=region)
    assert actual_hard_masked == expected_soft_masked
    assert actual_soft_masked == expected_soft_masked


def test_screen_pair_results(
    valid_primer_pairs: list[PrimerPair],
    genome_ref: Path,
    pair_primer_params: PrimerAndAmpliconParameters,
) -> None:
    """Test that `_has_acceptable_dinuc_run()` and `_screen_pair_results()` use
    `Primer3Parameters.primer_max_dinuc_bases` to disqualify primers when applicable.
    Create 2 sets of design input, the only difference being the length of allowable dinucleotide
    run in a primer (high_threshold = 6, low_threshold = 2).
    If one primer of a primer pair should have a dinucleotide run above the set threshold,
    then the pair is considered invalid."""
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
    design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=pair_primer_params,
        task=DesignPrimerPairsTask(),
    )

    lower_dinuc_thresh: PrimerAndAmpliconParameters = replace(
        pair_primer_params, primer_max_dinuc_bases=2
    )  # lower from 6 to 2
    altered_design_input: Primer3Input = Primer3Input(
        target=target,
        primer_and_amplicon_params=lower_dinuc_thresh,
        task=DesignPrimerPairsTask(),
    )
    with Primer3(genome_fasta=genome_ref) as designer:
        # all PrimerPairs have acceptable dinucleotide run lengths (threshold = 6)
        base_primer_pair_designs, base_dinuc_pair_failures = designer._screen_pair_results(
            design_input=design_input, designed_primer_pairs=valid_primer_pairs
        )
        assert len(base_dinuc_pair_failures) == 0
        assert design_input.primer_and_amplicon_params.primer_max_dinuc_bases is not None
        for primer_pair in base_primer_pair_designs:
            assert (
                primer_pair.left_primer.longest_dinucleotide_run_length()
                <= design_input.primer_and_amplicon_params.primer_max_dinuc_bases
            )
            assert (
                primer_pair.right_primer.longest_dinucleotide_run_length()
                <= design_input.primer_and_amplicon_params.primer_max_dinuc_bases
            )
            assert _has_acceptable_dinuc_run(
                design_input=design_input, oligo_design=primer_pair.left_primer
            )
            assert _has_acceptable_dinuc_run(
                design_input=design_input, oligo_design=primer_pair.right_primer
            )

        # 1 primer from every pair will fail lowered dinuc threshold of 2
        # As a result, no valid primer pairs will be emitted
        altered_designs, altered_dinuc_failures = designer._screen_pair_results(
            design_input=altered_design_input, designed_primer_pairs=valid_primer_pairs
        )
        assert altered_design_input.primer_and_amplicon_params is not None

        assert [
            design.longest_dinucleotide_run_length()
            > altered_design_input.primer_and_amplicon_params.primer_max_dinuc_bases
            for design in altered_dinuc_failures
        ]
        assert len(altered_designs) == 0


def test_build_failures(
    valid_primer_pairs: list[PrimerPair],
    genome_ref: Path,
    pair_primer_params: PrimerAndAmpliconParameters,
) -> None:
    """Test that `build_failures()` parses Primer3 `failure_strings` correctly and includes failures
    related to long dinucleotide runs."""
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)

    low_dinuc_thresh = replace(pair_primer_params, primer_max_dinuc_bases=2)  # lower from 6 to 2
    altered_design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=low_dinuc_thresh,
        task=DesignPrimerPairsTask(),
    )
    designer = Primer3(genome_fasta=genome_ref)
    primer_pair_designs, dinuc_pair_failures = designer._screen_pair_results(
        design_input=altered_design_input, designed_primer_pairs=valid_primer_pairs
    )
    # 3 primers fail for dinucleotide runs that are longer than `primer_max_dinuc_bases`
    assert len(dinuc_pair_failures) == 3
    test_failure_strings = [
        "considered 228, low tm 159, high tm 12, high hairpin stability 23, ok 34"
    ]
    fail_cases: list[Primer3Failure] = designer._build_failures(
        dinuc_pair_failures, test_failure_strings
    )
    # these fail_cases should match the cases in test_failure_strings excluding "considered" + "ok"
    for fail in fail_cases:
        if fail.reason == "low tm":
            assert fail.count == 159
        elif fail.reason == "high tm":
            assert fail.count == 12
        elif fail.reason == "high hairpin stability":
            assert fail.count == 23
        elif fail.reason == "long dinucleotide run":  # expect 3 failures due to dinucleotide run
            assert fail.count == 3


def test_build_failures_debugs(
    valid_primer_pairs: list[PrimerPair],
    genome_ref: Path,
    pair_primer_params: PrimerAndAmpliconParameters,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that we log a debug message in the event of an unknown Primer3Failure reason."""
    caplog.set_level(logging.DEBUG)
    target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)

    design_input = Primer3Input(
        target=target,
        primer_and_amplicon_params=pair_primer_params,
        task=DesignPrimerPairsTask(),
    )
    designer = Primer3(genome_fasta=genome_ref)
    primer_pair_designs, dinuc_pair_failures = designer._screen_pair_results(
        design_input=design_input, designed_primer_pairs=valid_primer_pairs
    )
    test_failure_strings = ["fabricated fail reason 1"]
    designer._build_failures(dinuc_pair_failures, test_failure_strings)
    expected_error_msg = "Unknown Primer3 failure reason"
    assert expected_error_msg in caplog.text


def test_primer3_result_primers_ok(
    valid_left_primers: list[Primer], valid_right_primers: list[Primer]
) -> None:
    primers: list[Primer] = valid_left_primers + valid_right_primers
    assert primers == Primer3Result(filtered_designs=primers, failures=[]).primers()


def test_primer3_result_primers_exception(valid_primer_pairs: list[PrimerPair]) -> None:
    result = Primer3Result(filtered_designs=valid_primer_pairs, failures=[])
    with pytest.raises(ValueError, match="Cannot call `primers` on `PrimerPair` results"):
        result.primers()


def test_primer3_result_as_primer_result_exception(valid_primer_pairs: list[PrimerPair]) -> None:
    result = Primer3Result(filtered_designs=valid_primer_pairs, failures=[])
    with pytest.raises(ValueError, match="Cannot call `as_primer_result` on `PrimerPair` results"):
        result.as_primer_result()


def test_primer3_result_primer_pairs_ok(valid_primer_pairs: list[PrimerPair]) -> None:
    assert valid_primer_pairs == (
        Primer3Result(filtered_designs=valid_primer_pairs, failures=[]).primer_pairs()
    )


def test_primer3_result_primer_pairs_exception(
    valid_left_primers: list[Primer], valid_right_primers: list[Primer]
) -> None:
    primers: list[Primer] = valid_left_primers + valid_right_primers
    result = Primer3Result(filtered_designs=primers, failures=[])
    with pytest.raises(ValueError, match="Cannot call `primer_pairs` on `Primer` results"):
        result.primer_pairs()


def test_primer3_result_as_primer_pair_result_exception(
    valid_left_primers: list[Primer], valid_right_primers: list[Primer]
) -> None:
    primers: list[Primer] = valid_left_primers + valid_right_primers
    result = Primer3Result(filtered_designs=primers, failures=[])
    with pytest.raises(ValueError, match="Cannot call `as_primer_pair_result` on `Primer` results"):
        result.as_primer_pair_result()


@pytest.mark.parametrize("max_amplicon_length", [100, 101])
def test_create_design_region(max_amplicon_length: int, genome_ref: Path) -> None:
    """If the target region is shorter than the max amplicon length, it should be padded to fit."""
    target_region = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)

    with Primer3(genome_fasta=genome_ref) as designer:
        design_region: Span = designer._create_design_region(
            target_region=target_region,
            max_amplicon_length=max_amplicon_length,
            min_primer_length=10,
        )

    assert design_region.length == 2 * max_amplicon_length - target_region.length


def test_create_design_region_raises_when_target_region_exceeds_max_amplicon_length(
    genome_ref: Path,
) -> None:
    """
    `_create_design_region()` should raise a ValueError when the target region is larger than the
    max amplicon length.
    """
    target_region = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)

    with Primer3(genome_fasta=genome_ref) as designer:
        with pytest.raises(ValueError, match="exceeds the maximum size"):
            designer._create_design_region(
                target_region=target_region, max_amplicon_length=10, min_primer_length=10
            )


def test_create_design_region_raises_when_primers_would_not_fit_in_design_region(
    genome_ref: Path,
) -> None:
    """
    `_create_design_region()` should raise a ValueError when the design region does not include
    sufficient space flanking the target for a primer to be designed. (i.e. when this space is less
    than the specified minimum primer length.)
    """
    target_region = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)

    with Primer3(genome_fasta=genome_ref) as designer:
        with pytest.raises(ValueError, match="exceeds the maximum size"):
            designer._create_design_region(
                target_region=target_region, max_amplicon_length=55, min_primer_length=10
            )
