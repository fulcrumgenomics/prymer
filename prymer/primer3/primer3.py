"""
# Primer3 Class and Methods

This module contains the [`Primer3`][prymer.primer3.primer3.Primer3] class, a class to
facilitate exchange of input and output data with Primer3, a command line tool.

Similar to the [`NtThermoAlign`][prymer.ntthal.NtThermoAlign] and
[`BwaAlnInteractive`][prymer.offtarget.bwa.BwaAlnInteractive] classes in the `prymer`
library, the Primer3 class extends the
[`ExecutableRunner`][prymer.util.executable_runner.ExecutableRunner] base class to
initiate an underlying subprocess, read and write input and output data, and gracefully terminate
any remaining subprocesses.

## Examples

The genome FASTA must be provided to the `Primer3` constructor, such that design and target
nucleotide sequences can be retrieved.  The full path to the `primer3` executable can provided,
otherwise it is assumed to be on the PATH.  Furthermore, optionally a
[`VariantLookup`][prymer.api.variant_lookup.VariantLookup] may be provided to
hard-mask the design and target regions as to avoid design primers over polymorphic sites.

```python
>>> from pathlib import Path
>>> from prymer.api.variant_lookup import VariantLookup, VariantOverlapDetector
>>> genome_fasta = Path("./tests/primer3/data/miniref.fa")
>>> genome_vcf = Path("./tests/primer3/data/miniref.variants.vcf.gz")
>>> variant_lookup: VariantLookup = VariantOverlapDetector(vcf_paths=[genome_vcf], min_maf=0.01, include_missing_mafs=False)
>>> designer = Primer3(genome_fasta=genome_fasta, variant_lookup=variant_lookup)

```

The `get_design_sequences()` method on `Primer3` is used to retrieve the soft and hard masked
sequences for a given region.  The hard-masked sequence replaces bases with `N`s that overlap
polymorphic sites found in the `VariantLookup` provided in the constructor.

```python
>>> design_region = Span(refname="chr2", start=9095, end=9120)
>>> soft_masked, hard_masked = designer.get_design_sequences(region=design_region)
>>> soft_masked
'AGTTACATTACAAAAGGCAGATTTCA'
>>> hard_masked
'AGTTANNNTACAAAAGGCAGATTTCA'

```

The `design_primers()` method on `Primer3` is used to design the primers given a
[`Primer3Input`][prymer.primer3.primer3_input.Primer3Input].  The latter includes all the
parameters and target region.

```python
>>> from prymer.primer3.primer3_parameters import Primer3Parameters
>>> from prymer.api import MinOptMax
>>> target = Span(refname="chr1", start=201, end=250, strand=Strand.POSITIVE)
>>> params = Primer3Parameters( \
    amplicon_sizes=MinOptMax(min=100, max=250, opt=200), \
    amplicon_tms=MinOptMax(min=55.0, max=100.0, opt=70.0), \
    primer_sizes=MinOptMax(min=29, max=31, opt=30), \
    primer_tms=MinOptMax(min=63.0, max=67.0, opt=65.0), \
    primer_gcs=MinOptMax(min=30.0, max=65.0, opt=45.0), \
)
>>> design_input = Primer3Input( \
    target=target, \
    params=params, \
    task=DesignLeftPrimersTask(), \
)
>>> left_result = designer.design_primers(design_input=design_input)

```

The `left_result` returns the [`Primer3Result`][prymer.primer3.primer3.Primer3Result]
container class.   It contains two attributes:
1. `filtered_designs`: filtered and ordered (by objective function score) list of primer pairs or
    single primers that were returned by Primer3.
2. `failures`: ordered list of [`Primer3Failures`][prymer.primer3.primer3.Primer3Failure]
    detailing design failure reasons and corresponding count.

In this case, there are two failures reasons:

```python
>>> for failure in left_result.failures: \
    print(failure)
Primer3Failure(reason=<Primer3FailureReason.HIGH_TM: 'high tm'>, count=406)
Primer3Failure(reason=<Primer3FailureReason.GC_CONTENT: 'GC content failed'>, count=91)

```

While`filtered_designs` attribute on `Primer3Result` may be used to access the list of primers or
primer pairs, it is more convenient to use the `primers()` and `primer_pairs()` methods
to return the designed primers or primer pairs (use the method corresponding to the input task) so
that the proper type is returned (i.e. [`Primer`][prymer.api.primer.Primer] or
[`PrimerPair`][prymer.api.primer_pair.PrimerPair]).

```python
>>> for primer in left_result.primers(): \
    print(primer)
TCTGAACAGGACGAACTGGATTTCCTCAT	65.686	1.953897	chr1:163-191:+
CTCTGAACAGGACGAACTGGATTTCCTCAT	66.152	2.293213	chr1:162-191:+
TCTGAACAGGACGAACTGGATTTCCTCATG	66.33	2.514048	chr1:163-192:+
AACAGGACGAACTGGATTTCCTCATGGAA	66.099	2.524986	chr1:167-195:+
CTGAACAGGACGAACTGGATTTCCTCATG	65.47	2.556859	chr1:164-192:+

```

Finally, the designer should be closed to terminate the sub-process:

```python
>>> designer.close()
True

```

`Primer3` is also context a manager, and so can be used with a `with` clause:

```python
>>> with Primer3(genome_fasta=genome_fasta) as designer: \
    pass  # use designer here!

```

"""  # noqa: E501

import logging
import subprocess
import typing
from collections import Counter
from dataclasses import dataclass
from dataclasses import replace
from pathlib import Path
from typing import Generic
from typing import Optional
from typing import TypeVar
from typing import Union
from typing import assert_never

import pysam
from fgpyo import sam
from fgpyo.fasta.sequence_dictionary import SequenceDictionary
from fgpyo.sam import reader
from fgpyo.sequence import reverse_complement
from fgpyo.util.metric import Metric

from prymer.api.primer import Primer
from prymer.api.primer_like import PrimerLike
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.api.span import Strand
from prymer.api.variant_lookup import SimpleVariant
from prymer.api.variant_lookup import VariantLookup
from prymer.primer3.primer3_failure_reason import Primer3FailureReason
from prymer.primer3.primer3_input import Primer3Input
from prymer.primer3.primer3_input_tag import Primer3InputTag
from prymer.primer3.primer3_task import DesignLeftPrimersTask
from prymer.primer3.primer3_task import DesignPrimerPairsTask
from prymer.primer3.primer3_task import DesignRightPrimersTask
from prymer.util.executable_runner import ExecutableRunner


@dataclass(init=True, slots=True, frozen=True)
class Primer3Failure(Metric["Primer3Failure"]):
    """Encapsulates how many designs failed for a given reason.
    Extends the `fgpyo.util.metric.Metric` class, which will facilitate writing out results for
    primer design QC etc.

    Attributes:
        reason: the reason the design failed
        count: how many designs failed
    """

    reason: Primer3FailureReason
    count: int


PrimerLikeType = TypeVar("PrimerLikeType", bound=PrimerLike)
"""Type variable for a `Primer3Result`, which must implement `PrimerLike`"""


@dataclass(init=True, slots=True, frozen=True)
class Primer3Result(Generic[PrimerLikeType]):
    """Encapsulates Primer3 design results (both valid designs and failures).

    Attributes:
        filtered_designs: filtered and ordered (by objective function score) list of primer
            pairs or single primers that were returned by Primer3
        failures: ordered list of Primer3Failures detailing design failure reasons and corresponding
            count
    """

    filtered_designs: list[PrimerLikeType]
    failures: list[Primer3Failure]

    def as_primer_result(self) -> "Primer3Result[Primer]":
        """Returns this Primer3Result assuming the design results are of type `Primer`."""
        if len(self.filtered_designs) > 0 and not isinstance(self.filtered_designs[0], Primer):
            raise ValueError("Cannot call `as_primer_result` on `PrimerPair` results")
        return typing.cast(Primer3Result[Primer], self)

    def as_primer_pair_result(self) -> "Primer3Result[PrimerPair]":
        """Returns this Primer3Result assuming the design results are of type `PrimerPair`."""
        if len(self.filtered_designs) > 0 and not isinstance(self.filtered_designs[0], PrimerPair):
            raise ValueError("Cannot call `as_primer_pair_result` on `Primer` results")
        return typing.cast(Primer3Result[PrimerPair], self)

    def primers(self) -> list[Primer]:
        """Returns the design results as a list `Primer`s"""
        try:
            return self.as_primer_result().filtered_designs
        except ValueError as ex:
            raise ValueError("Cannot call `primers` on `PrimerPair` results") from ex

    def primer_pairs(self) -> list[PrimerPair]:
        """Returns the design results as a list `PrimerPair`s"""
        try:
            return self.as_primer_pair_result().filtered_designs
        except ValueError as ex:
            raise ValueError("Cannot call `primer_pairs` on `Primer` results") from ex


class Primer3(ExecutableRunner):
    """
    Enables interaction with command line tool, primer3.

    Attributes:
        _fasta: file handle to the open reference genome file
        _dict: the sequence dictionary that corresponds to the provided reference genome file
    """

    def __init__(
        self,
        genome_fasta: Path,
        executable: Optional[str] = None,
        variant_lookup: Optional[VariantLookup] = None,
    ) -> None:
        """
        Args:
            genome_fasta: Path to reference genome .fasta file
            executable: string representation of the path to primer3_core
            variant_lookup: VariantLookup object to facilitate hard-masking variants

        Assumes the sequence dictionary is located adjacent to the .fasta file and has the same
        base name with a .dict suffix.

        """
        executable_path = ExecutableRunner.validate_executable_path(
            executable="primer3_core" if executable is None else executable
        )
        command: list[str] = [f"{executable_path}"]

        self.variant_lookup = variant_lookup
        self._fasta = pysam.FastaFile(filename=f"{genome_fasta}")

        dict_path = genome_fasta.with_suffix(".dict")
        # TODO: This is a placeholder while waiting for #160  to be resolved
        # https://github.com/fulcrumgenomics/fgpyo/pull/160
        with reader(dict_path, file_type=sam.SamFileType.SAM) as fh:
            self._dict: SequenceDictionary = SequenceDictionary.from_sam(header=fh.header)

        super().__init__(command=command, stderr=subprocess.STDOUT)

    def close(self) -> bool:
        """Closes fasta file regardless of underlying subprocess status.
        Logs an error if the underlying subprocess is not successfully closed.

        Returns:
            True: if the subprocess was terminated successfully
            False: if the subprocess failed to terminate or was not already running
        """
        self._fasta.close()
        subprocess_close = super().close()
        if not subprocess_close:
            logging.debug("Did not successfully close underlying subprocess")
        return subprocess_close

    def get_design_sequences(self, region: Span) -> tuple[str, str]:
        """Extracts the reference sequence that corresponds to the design region.

        Args:
            region: the region of the genome to be extracted

        Returns:
            A tuple of two sequences: the sequence for the region, and the sequence for the region
            with variants hard-masked as Ns

        """
        # pysam.fetch: 0-based, half-open intervals
        soft_masked = self._fasta.fetch(
            reference=region.refname, start=region.start - 1, end=region.end
        )

        if self.variant_lookup is None:
            hard_masked = soft_masked
            return soft_masked, hard_masked

        overlapping_variants: list[SimpleVariant] = self.variant_lookup.query(
            refname=region.refname, start=region.start, end=region.end
        )
        positions: list[int] = []
        for variant in overlapping_variants:
            # FIXME
            positions.extend(range(variant.pos, variant.end + 1))

        filtered_positions = [pos for pos in positions if region.start <= pos <= region.end]
        soft_masked_list = list(soft_masked)
        for pos in filtered_positions:
            soft_masked_list[region.get_offset(pos)] = (
                "N"  # get relative coord of filtered position and mask to N
            )
        # convert list back to string
        hard_masked = "".join(soft_masked_list)
        return soft_masked, hard_masked

    @staticmethod
    def _is_valid_primer(design_input: Primer3Input, primer_design: Primer) -> bool:
        return (
            primer_design.longest_dinucleotide_run_length()
            <= design_input.params.primer_max_dinuc_bases
        )

    @staticmethod
    def _screen_pair_results(
        design_input: Primer3Input, designed_primer_pairs: list[PrimerPair]
    ) -> tuple[list[PrimerPair], list[Primer]]:
        """Screens primer pair designs emitted by Primer3 for dinucleotide run length.

        Args:
            design_input: the target region, design task, specifications, and scoring penalties
            designed_primer_pairs: the unfiltered primer pair designs emitted by Primer3

        Returns:
            valid_primer_pair_designs: primer pairs within specifications
            dinuc_pair_failures: single primer designs that failed the `max_dinuc_bases` threshold
        """
        valid_primer_pair_designs: list[PrimerPair] = []
        dinuc_pair_failures: list[Primer] = []
        for primer_pair in designed_primer_pairs:
            valid: bool = True
            if (
                primer_pair.left_primer.longest_dinucleotide_run_length()
                > design_input.params.primer_max_dinuc_bases
            ):  # if the left primer has too many dinucleotide bases, fail it
                dinuc_pair_failures.append(primer_pair.left_primer)
                valid = False
            if (
                primer_pair.right_primer.longest_dinucleotide_run_length()
                > design_input.params.primer_max_dinuc_bases
            ):  # if the right primer has too many dinucleotide bases, fail it
                dinuc_pair_failures.append(primer_pair.right_primer)
                valid = False
            if valid:  # if neither failed, append the pair to a list of valid designs
                valid_primer_pair_designs.append(primer_pair)
        return valid_primer_pair_designs, dinuc_pair_failures

    def design_primers(self, design_input: Primer3Input) -> Primer3Result:  # noqa: C901
        """Designs primers or primer pairs given a target region.

        Args:
            design_input: encapsulates the target region, design task, specifications, and scoring
                penalties

        Returns:
            Primer3Result containing both the valid and failed designs emitted by Primer3

        Raises:
            RuntimeError: if underlying subprocess is not alive
            ValueError: if Primer3 returns errors or does not return output
            ValueError: if Primer3 output is malformed
            ValueError: if an unknown design task is given
        """

        if not self.is_alive:
            raise RuntimeError(
                f"Error, trying to use a subprocess that has already been "
                f"terminated, return code {self._subprocess.returncode}"
            )

        design_region: Span = self._create_design_region(
            target_region=design_input.target,
            max_amplicon_length=design_input.params.max_amplicon_length,
            min_primer_length=design_input.params.min_primer_length,
        )

        soft_masked, hard_masked = self.get_design_sequences(design_region)
        global_primer3_params = {
            Primer3InputTag.PRIMER_FIRST_BASE_INDEX: 1,
            Primer3InputTag.PRIMER_EXPLAIN_FLAG: 1,
            Primer3InputTag.SEQUENCE_TEMPLATE: hard_masked,
        }

        assembled_primer3_tags = {
            **global_primer3_params,
            **design_input.to_input_tags(design_region=design_region),
        }

        # Submit inputs to primer3
        for tag, value in assembled_primer3_tags.items():
            self._subprocess.stdin.write(f"{tag}={value}")
            self._subprocess.stdin.write("\n")
        self._subprocess.stdin.write("=\n")
        self._subprocess.stdin.flush()

        error_lines: list[str] = []  # list of errors as reported by primer3
        primer3_results: dict[str, str] = {}  # key-value pairs of results reported by Primer3

        def primer3_error(message: str) -> None:
            """Formats the Primer3 error and raises a ValueError."""
            error_message = f"{message}: "
            # add in any reported PRIMER_ERROR
            if "PRIMER_ERROR" in primer3_results:
                error_message += primer3_results["PRIMER_ERROR"]
            # add in any error lines
            if len(error_lines) > 0:
                error_message += "\n".join(f"\t\t{e}" for e in error_lines)
            # raise the exception now
            raise ValueError(error_message)

        while True:
            # Get the next line.  Since we want to distinguish between empty lines, which we ignore,
            # and the end-of-file, which is just an empty string, check for an empty string before
            # stripping the line of any trailing newline or carriage return characters.
            line: str = self._subprocess.stdout.readline()
            if line == "":  # EOF
                primer3_error("Primer3 exited prematurely")
            line = line.rstrip("\r\n")

            if line == "=":  # stop when we find the line just "="
                break
            elif line == "":  # ignore empty lines
                continue
            elif "=" not in line:  # error lines do not have the equals character in them, usually
                error_lines.append(line)
            else:  # parse and store the result
                key, value = line.split("=", maxsplit=1)
                # Because Primer3 will emit both the input given and the output generated, we
                # discard the input that is echo'ed back by looking for tags (keys)
                # that do not match any Primer3InputTag
                if not any(key == item.value for item in Primer3InputTag):
                    primer3_results[key] = value

        # Check for any errors.  Typically, these are in error_lines, but also the results can
        # contain the PRIMER_ERROR key.
        if "PRIMER_ERROR" in primer3_results or len(error_lines) > 0:
            primer3_error("Primer3 failed")

        match design_input.task:
            case DesignPrimerPairsTask():  # Primer pair design
                all_pair_results: list[PrimerPair] = Primer3._build_primer_pairs(
                    design_input=design_input,
                    design_results=primer3_results,
                    design_region=design_region,
                    unmasked_design_seq=soft_masked,
                )
                return Primer3._assemble_primer_pairs(
                    design_input=design_input,
                    design_results=primer3_results,
                    unfiltered_designs=all_pair_results,
                )

            case DesignLeftPrimersTask() | DesignRightPrimersTask():  # Single primer design
                all_single_results = Primer3._build_primers(
                    design_input=design_input,
                    design_results=primer3_results,
                    design_region=design_region,
                    design_task=design_input.task,
                    unmasked_design_seq=soft_masked,
                )
                return Primer3._assemble_primers(
                    design_input=design_input,
                    design_results=primer3_results,
                    unfiltered_designs=all_single_results,
                )

            case _ as unreachable:
                assert_never(unreachable)

    @staticmethod
    def _build_primers(
        design_input: Primer3Input,
        design_results: dict[str, str],
        design_region: Span,
        design_task: Union[DesignLeftPrimersTask, DesignRightPrimersTask],
        unmasked_design_seq: str,
    ) -> list[Primer]:
        """
        Builds a list of left or right primers from Primer3 output.

        Args:
            design_input: the target region, design task, specifications, and scoring penalties
            design_results: design results emitted by Primer3 and captured by design_primers()
            design_region: the padded design region
            design_task: the design task
            unmasked_design_seq: the reference sequence corresponding to the target region

        Returns:
            primers: a list of unsorted and unfiltered primer designs emitted by Primer3

        Raises:
            ValueError: if Primer3 does not return primer designs
        """
        count_tag = design_input.task.count_tag

        maybe_count: Optional[str] = design_results.get(count_tag)
        if maybe_count is None:  # no count tag was found
            if "PRIMER_ERROR" in design_results:
                primer_error = design_results["PRIMER_ERROR"]
                raise ValueError(f"Primer3 returned an error: {primer_error}")
            else:
                raise ValueError(f"Primer3 did not return the count tag: {count_tag}")
        count: int = int(maybe_count)

        primers = []
        for idx in range(count):
            key = f"PRIMER_{design_task.task_type}_{idx}"
            str_position, str_length = design_results[key].split(",", maxsplit=1)
            position, length = int(str_position), int(str_length)  # position is 1-based

            match design_task:
                case DesignLeftPrimersTask():
                    span = design_region.get_subspan(
                        offset=position - 1, subspan_length=length, strand=Strand.POSITIVE
                    )
                case DesignRightPrimersTask():
                    start = position - length + 1  # start is 1-based
                    span = design_region.get_subspan(
                        offset=start - 1, subspan_length=length, strand=Strand.NEGATIVE
                    )
                case _ as unreachable:
                    assert_never(unreachable)  # pragma: no cover

            slice_offset = design_region.get_offset(span.start)
            slice_end = design_region.get_offset(span.end) + 1

            # remake the primer sequence from the un-masked genome sequence just in case
            bases = unmasked_design_seq[slice_offset:slice_end]
            if span.strand == Strand.NEGATIVE:
                bases = reverse_complement(bases)

            primers.append(
                Primer(
                    bases=bases,
                    tm=float(design_results[f"PRIMER_{design_task.task_type}_{idx}_TM"]),
                    penalty=float(design_results[f"PRIMER_{design_task.task_type}_{idx}_PENALTY"]),
                    span=span,
                )
            )
        return primers

    @staticmethod
    def _assemble_primers(
        design_input: Primer3Input, design_results: dict[str, str], unfiltered_designs: list[Primer]
    ) -> Primer3Result:
        """Helper function to organize primer designs into valid and failed designs.

        Wraps `Primer3._is_valid_primer()` and `Primer3._build_failures()` to filter out designs
        with dinucleotide runs that are too long and extract additional failure reasons emitted by
        Primer3.

        Args:
            design_input: encapsulates the target region, design task, specifications,
             and scoring penalties
            unfiltered_designs: list of primers emitted from Primer3
             design_results: key-value pairs of results reported by Primer3

        Returns:
            primer_designs: a `Primer3Result` that encapsulates valid and failed designs
        """
        valid_primer_designs = [
            design
            for design in unfiltered_designs
            if Primer3._is_valid_primer(primer_design=design, design_input=design_input)
        ]
        dinuc_failures = [
            design
            for design in unfiltered_designs
            if not Primer3._is_valid_primer(primer_design=design, design_input=design_input)
        ]

        failure_strings = [design_results[f"PRIMER_{design_input.task.task_type}_EXPLAIN"]]
        failures = Primer3._build_failures(dinuc_failures, failure_strings)
        primer_designs: Primer3Result = Primer3Result(
            filtered_designs=valid_primer_designs, failures=failures
        )
        return primer_designs

    @staticmethod
    def _build_primer_pairs(
        design_input: Primer3Input,
        design_results: dict[str, str],
        design_region: Span,
        unmasked_design_seq: str,
    ) -> list[PrimerPair]:
        """
        Builds a list of primer pairs from single primer designs emitted from Primer3.

        Args:
            design_input: the target region, design task, specifications, and scoring penalties
            design_results: design results emitted by Primer3 and captured by design_primers()
            design_region: the padded design region
            unmasked_design_seq: the reference sequence corresponding to the target region

        Returns:
            primer_pairs: a list of unsorted and unfiltered paired primer designs emitted by Primer3

        Raises:
            ValueError: if Primer3 does not return the same number of left and right designs
        """
        left_primers = Primer3._build_primers(
            design_input=design_input,
            design_results=design_results,
            design_region=design_region,
            design_task=DesignLeftPrimersTask(),
            unmasked_design_seq=unmasked_design_seq,
        )

        right_primers = Primer3._build_primers(
            design_input=design_input,
            design_results=design_results,
            design_region=design_region,
            design_task=DesignRightPrimersTask(),
            unmasked_design_seq=unmasked_design_seq,
        )

        def _build_primer_pair(num: int, primer_pair: tuple[Primer, Primer]) -> PrimerPair:
            """Builds the `PrimerPair` object from input left and right primers."""
            left_primer = primer_pair[0]
            right_primer = primer_pair[1]
            amplicon = replace(left_primer.span, end=right_primer.span.end)
            slice_offset = design_region.get_offset(amplicon.start)
            slice_end = slice_offset + amplicon.length

            return PrimerPair(
                left_primer=left_primer,
                right_primer=right_primer,
                amplicon_tm=float(design_results[f"PRIMER_PAIR_{num}_PRODUCT_TM"]),
                penalty=float(design_results[f"PRIMER_PAIR_{num}_PENALTY"]),
                amplicon_sequence=unmasked_design_seq[slice_offset:slice_end],
            )

        #  Primer3 returns an equal number of left and right primers during primer pair design
        if len(left_primers) != len(right_primers):
            raise ValueError("Primer3 returned a different number of left and right primers.")
        primer_pairs: list[PrimerPair] = [
            _build_primer_pair(num, primer_pair)
            for num, primer_pair in enumerate(zip(left_primers, right_primers, strict=True))
        ]
        return primer_pairs

    @staticmethod
    def _assemble_primer_pairs(
        design_input: Primer3Input,
        design_results: dict[str, str],
        unfiltered_designs: list[PrimerPair],
    ) -> Primer3Result:
        """Helper function to organize primer pairs into valid and failed designs.

        Wraps `Primer3._screen_pair_results()` and `Primer3._build_failures()` to filter out designs
        with dinucleotide runs that are too long and extract additional failure reasons emitted by
        Primer3.

        Args:
            design_input: encapsulates the target region, design task, specifications,
             and scoring penalties
            unfiltered_designs: list of primer pairs emitted from Primer3
             design_results: key-value pairs of results reported by Primer3

        Returns:
            primer_designs: a `Primer3Result` that encapsulates valid and failed designs
        """
        valid_primer_pair_designs: list[PrimerPair]
        dinuc_pair_failures: list[Primer]
        valid_primer_pair_designs, dinuc_pair_failures = Primer3._screen_pair_results(
            design_input=design_input, designed_primer_pairs=unfiltered_designs
        )

        failure_strings = [
            design_results["PRIMER_PAIR_EXPLAIN"],
            design_results["PRIMER_LEFT_EXPLAIN"],
            design_results["PRIMER_RIGHT_EXPLAIN"],
        ]
        pair_failures = Primer3._build_failures(dinuc_pair_failures, failure_strings)
        primer_designs = Primer3Result(
            filtered_designs=valid_primer_pair_designs, failures=pair_failures
        )

        return primer_designs

    @staticmethod
    def _build_failures(
        dinuc_failures: list[Primer],
        failure_strings: list[str],
    ) -> list[Primer3Failure]:
        """Extracts the reasons why designs that were considered by Primer3 failed
         (when there were failures).

        The set of failures is returned sorted from those with most
        failures to those with least.

        Args:
            dinuc_failures: primer designs with a dinucleotide run longer than the allowed maximum
            failure_strings: explanations (strings) emitted by Primer3 about failed designs


        Returns:
            a list of Primer3Failure objects
        """

        by_fail_count: Counter[Primer3FailureReason] = Primer3FailureReason.parse_failures(
            *failure_strings
        )
        # Count how many individual primers failed for dinuc runs
        num_dinuc_failures = len(set(dinuc_failures))
        if num_dinuc_failures > 0:
            by_fail_count[Primer3FailureReason.LONG_DINUC] = num_dinuc_failures
        return [Primer3Failure(reason, count) for reason, count in by_fail_count.most_common()]

    def _create_design_region(
        self,
        target_region: Span,
        max_amplicon_length: int,
        min_primer_length: int,
    ) -> Span:
        """
        Construct a design region surrounding the target region.

        The target region is padded on both sides by the maximum amplicon length, minus the length
        of the target region itself.

        If the target region cannot be padded by at least the minimum primer length on both sides,
        a `ValueError` is raised.

        Raises:
            ValueError: If the target region is too large to be padded.

        """
        # Pad the target region on both sides by the maximum amplicon length (minus the length of
        # the target). This ensures that the design region covers the complete window of potentially
        # valid primer pairs.
        padding: int = max_amplicon_length - target_region.length

        # Apply the padding, ensuring that we don't run out-of-bounds on the target contig.
        contig_length: int = self._dict[target_region.refname].length
        design_start: int = max(1, target_region.start - padding)
        design_end: int = min(target_region.end + padding, contig_length)

        # Validate that our design window includes sufficient space for a primer to be designed on
        # each side of the target region.
        left_design_window: int = target_region.start - design_start
        right_design_window: int = design_end - target_region.end
        if left_design_window < min_primer_length or right_design_window < min_primer_length:
            raise ValueError(
                f"Target region {target_region} exceeds the maximum size compatible with a "
                f"maximum amplicon length of {max_amplicon_length} and a minimum primer length of "
                f"{min_primer_length}. The maximum amplicon length should exceed the length of "
                "the target region by at least twice the minimum primer length."
            )

        # Return the validated design region.
        design_region: Span = replace(
            target_region,
            start=design_start,
            end=design_end,
        )

        return design_region
