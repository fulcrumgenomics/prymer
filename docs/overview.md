# Overview

The `prymer` Python library is intended to be used for three main purposes:

1. [Clustering targets](#clustering-targets) into larger amplicons prior to designing primers.
2. [Designing primers](#designing-primers) (left or right) or primer pairs using Primer3 for each target from (1).
3. [Build and Picking a set of primer pairs](#build-and-picking-primer-pairs) from the designed primer pairs produced in (2).

## Clustering Targets

Optionally, input targets may be clustered into larger amplicons prior to designing primers.  These amplicons wholly 
contain the input targets, and are used as input for primer design.  The [`cluster_intervals()`][prymer.api.clustering.cluster_intervals] 
method in the [`prymer.api.clustering`][prymer.api.clustering] module is used to cluster targets into larger 
amplicons prior to primer design.

## Designing Primers

Designing primers (left or right) or primer pairs using Primer3 is primarily performed using the 
[`Primer3`][prymer.primer3.primer3.Primer3] class, which wraps the
[`primer3` command line tool](https://github.com/primer3-org/primer3).  The 
[`design_oligos()`][prymer.primer3.primer3.Primer3.design_oligos] facilitates the design of primers (single and paired) and/or internal probes 
for a single target. The `Primer3` instance is intended to be re-used to design primers across multiple targets, or 
re-design (after changing parameters) for the same target, or both!

Common input parameters for designing primers are specified in [`Primer3Parameters()`][prymer.primer3.primer3_parameters.Primer3Parameters] and 
[`PrimerAndAmpliconWeights()`][prymer.primer3.primer3_weights.PrimerAndAmpliconWeights], while the task type (left primer,
right primer, or primer pair design) is specified with the corresponding 
[`Primer3Task`][prymer.primer3.primer3_task.Primer3Task]. Penalty weights for designing internal probes are specified in [`ProbeWeights()`][prymer.primer3.primer3_weights.ProbeWeights]

The result of a primer design is encapsulated in the [`Primer3Result`][prymer.primer3.primer3.Primer3Result] class.  It
provides the primers (or primer pairs) that were designed, as well as a list of reasons some primers were not returned, 
for example exceeding the melting temperature threshold, too high GC content, and so on.  These failures are 
encapsulated in the [`Primer3Failures`][prymer.primer3.primer3.Primer3Failure] class.

The [`Primer3Result`][prymer.primer3.primer3.Primer3Result] returned by the primer design contains either a list of 
[`Primer`][prymer.api.primer.Primer]s or [`PrimerPair`][prymer.api.primer_pair.PrimerPair]s, depending on the 
[`Primer3Task`][prymer.primer3.primer3_task.Primer3Task] specified in the input parameters.
These can be subsequently filtered or examined.

## Build and Picking Primer Pairs

It is recommended to design left and right primers individually, then to screen each primer for off-target mappings
using the [`OffTargetDetector`][prymer.offtarget.offtarget_detector.OffTargetDetector] class, which wraps 
`bwa aln` to identify off-target primer mappings.  

The remaining primers may then be used with the 
[`build_primer_pairs()`][prymer.api.picking.build_primer_pairs] method to build primer pairs
from all combinations of the left and right primers.
The produced primer pairs are scored in a manner similar to Primer3 using the [`score()`][prymer.api.picking.score] method.
The [`FilterParams`][prymer.api.picking.FilteringParams] class is used to provide parameters for scoring.

Next, the [`pick_top_primer_pairs()`][prymer.api.picking.pick_top_primer_pairs] method is used to select up to
a maximum number of primer pairs.  The primer pairs are selected in the order of lowest penalty (highest score).  As
part of this process, each primer pair must:

1. Have an amplicon in the desired size range (see [`is_acceptable_primer_pair`][prymer.api.picking.is_acceptable_primer_pair]).
2. Have an amplicon melting temperature in the desired range (see [`is_acceptable_primer_pair`][prymer.api.picking.is_acceptable_primer_pair]).
3. Not have too many off-targets (see [`OffTargetDetector.check_one()`][prymer.offtarget.offtarget_detector.OffTargetDetector.check_one]).
4. Not have primer pairs that overlap too much (see [`check_primer_overlap()`][prymer.api.picking.check_primer_overlap]).
5. Not form a dimer with a melting temperature above a specified threshold (see the [`max_dimer_tm` attribute in `FilterParams`][prymer.api.picking.FilteringParams]).

Checking for dimers may be performed using the [`NtThermoAlign`][prymer.ntthal.NtThermoAlign] command line executable,
and can be passed to [`pick_top_primer_pairs()`][prymer.api.picking.pick_top_primer_pairs] as follows:

```python
from prymer.ntthal import NtThermoAlign
from prymer.api.picking import FilteringParams, pick_top_primer_pairs
params = FilteringParams(...)
dimer_checker = NtThermoAlign()
pick_top_primer_pairs(
    is_dimer_tm_ok=lambda s1, s2: (
        dimer_checker.duplex_tm(s1=s1, s2=s2) <= params.max_dimer_tm
    ),
    ...
)

```

For convenience, the [`build_and_pick_primer_pairs()`][prymer.api.picking.build_and_pick_primer_pairs] method combines
both the [`build_primer_pairs()`][prymer.api.picking.build_primer_pairs] and 
[`pick_top_primer_pairs()`][prymer.api.picking.pick_top_primer_pairs] methods in single invocation.

The resulting primer pairs may be further examined.