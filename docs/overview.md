# Overview

The `prymer` Python library is intended to be used for three main purposes:

1. [Clustering targets](#clustering-targets) into larger amplicons prior to designing primers.
2. [Designing](#designing-primers) primers (single and paired) and internal hybridization probes using Primer3 for each target from (1).
3. [Build and Picking a set of primer pairs](#build-and-picking-primer-pairs) from the design candidates produced in (2).

## Clustering Targets

Optionally, input targets may be clustered into larger amplicons prior to designing primers.  These amplicons wholly 
contain the input targets, and are used as input for primer design.  The [`cluster_intervals()`][prymer.api.clustering.cluster_intervals] 
method in the [`prymer.api.clustering`][prymer.api.clustering] module is used to cluster targets into larger 
amplicons prior to primer design.

## Designing Primers

Designing primers (left or right) or primer pairs using Primer3 is primarily performed using the 
[`Primer3`][prymer.primer3.primer3.Primer3] class, which wraps the
[`primer3` command line tool](https://github.com/primer3-org/primer3).  The 
[`design()`][prymer.primer3.primer3.Primer3.design] method facilitates the design of primers (single and paired) and internal hybridization probes 
for a single target. The `Primer3` instance is intended to be re-used to design primers across multiple targets, or 
re-design (after changing parameters) for the same target, or both!

Common input parameters for designing primers are specified in [`PrimerAndAmpliconParameters()`][prymer.primer3.primer3_parameters.PrimerAndAmpliconParameters] and 
[`PrimerAndAmpliconWeights()`][prymer.primer3.primer3_weights.PrimerAndAmpliconWeights], while the task type (left primer,
right primer, or primer pair design) is specified with the corresponding 
[`Primer3Task`][prymer.primer3.primer3_task.Primer3Task]. 
Design specifications for designing probes are stored in [`ProbeParameters()`][prymer.primer3.primer3_parameters.ProbeParameters]. 
Penalty weights for designing internal probes are specified in [`ProbeWeights()`][prymer.primer3.primer3_weights.ProbeWeights].

The result of a primer design is encapsulated in the [`Primer3Result`][prymer.primer3.primer3.Primer3Result] class.  It
provides the primers, probes, or primer pairs that were designed, as well as a list of reasons some primers were not returned, 
for example exceeding the melting temperature threshold, too high GC content, and so on.  These failures are 
encapsulated in the [`Primer3Failures`][prymer.primer3.primer3.Primer3Failure] class.

The [`Primer3Result`][prymer.primer3.primer3.Primer3Result] returned by the primer design contains either a list of 
[`Oligo`][prymer.api.primer.Oligo]s or [`PrimerPair`][prymer.api.primer_pair.PrimerPair]s, depending on the 
[`Primer3Task`][prymer.primer3.primer3_task.Primer3Task] specified in the input parameters.
These can be subsequently filtered or examined.

## Build and Picking Primer Pairs

It is recommended to design left and right primers individually, then to screen each primer for off-target mappings
using the [`OffTargetDetector`][prymer.offtarget.offtarget_detector.OffTargetDetector] class, which wraps 
`bwa aln` to identify off-target primer mappings.  

The remaining primers may then be used with the 
[`build_primer_pairs()`][prymer.api.picking.build_primer_pairs] method to build primer pairs
from all combinations of the left and right primers that are within acceptable amplicon size and tm ranges.
Optionally, if a `max_heterodimer_tm` is provided, primer pairs will be screened to ensure that the left and right
primers do not dimerize with a Tm greater than the one provided.
The produced primer pairs are scored in a manner similar to Primer3 using the [`score()`][prymer.api.picking.score] method.