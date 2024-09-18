"""
# Methods and classes to detect off-target mappings for primers and primer pairs

The [`OffTargetDetector`][prymer.offtarget.offtarget_detector.OffTargetDetector] uses
[`BwaAlnInteractive`][prymer.offtarget.bwa.BwaAlnInteractive] to search for off targets for
one or more primers or primer pairs.

The [`filter()`][prymer.offtarget.offtarget_detector.OffTargetDetector.filter] method
filters an iterable of `Primers` to return only those that have less than a given maximum number of off-target hits
to the genome.

```python
>>> from pathlib import Path
>>> from prymer.api.span import Strand
>>> ref_fasta = Path("./tests/offtarget/data/miniref.fa")
>>> left_primer =  Primer(bases="AAAAA", tm=37, penalty=0, span=Span("chr1", start=67, end=71))
>>> right_primer = Primer(bases="TTTTT", tm=37, penalty=0, span=Span("chr1", start=75, end=79, strand=Strand.NEGATIVE))
>>> detector = OffTargetDetector(ref=ref_fasta, max_primer_hits=204, max_primer_pair_hits=1, three_prime_region_length=20, max_mismatches_in_three_prime_region=0, max_mismatches=0, max_amplicon_size=250)
>>> len(detector.filter(primers=[left_primer, right_primer]))  # keep all
2
>>> detector.close()
>>> detector = OffTargetDetector(ref=ref_fasta, max_primer_hits=203, max_primer_pair_hits=1, three_prime_region_length=20, max_mismatches_in_three_prime_region=0, max_mismatches=0, max_amplicon_size=250)
>>> len(detector.filter(primers=[left_primer, right_primer]))  # keep  none
0
>>> detector.close()

```

The [`check_one()`][prymer.offtarget.offtarget_detector.OffTargetDetector.check_one] and
the [`check_all()`][prymer.offtarget.offtarget_detector.OffTargetDetector.check_all]
methods check one or multiple primer pairs respectively for off-target mappings, returning an
[`OffTargetResult`][prymer.offtarget.offtarget_detector.OffTargetResult], or mapping from
each `PrimerPair` to its corresponding `OffTargetResult`.  This result contains information about
the off-target mappings.

```python
>>> detector = OffTargetDetector(ref=ref_fasta, max_primer_hits=1, max_primer_pair_hits=1, three_prime_region_length=20, max_mismatches_in_three_prime_region=0, max_mismatches=0, max_amplicon_size=250)
>>> primer_pair = PrimerPair(left_primer=left_primer, right_primer=right_primer, amplicon_tm=37, penalty=0.0)
>>> detector.check_one(primer_pair=primer_pair)
OffTargetResult(primer_pair=..., passes=False, cached=False, spans=[], left_primer_spans=[], right_primer_spans=[])
>>> list(detector.check_all(primer_pairs=[primer_pair]).values())
[OffTargetResult(primer_pair=..., passes=False, cached=True, spans=[], left_primer_spans=[], right_primer_spans=[])]
>>> detector = OffTargetDetector(ref=ref_fasta, max_primer_hits=204, max_primer_pair_hits=856, three_prime_region_length=20, max_mismatches_in_three_prime_region=0, max_mismatches=0, max_amplicon_size=250)
>>> result = detector.check_one(primer_pair=primer_pair)
>>> len(result.spans)
856
>>> len(result.left_primer_spans)
204
>>> len(result.right_primer_spans)
204

```

Finally, the [`mappings_of()`][prymer.offtarget.offtarget_detector.OffTargetDetector.mappings_of]
method maps individual primers (`Primer`s).

```python
>>> p1: Primer = Primer(tm=37, penalty=0, span=Span(refname="chr1", start=1, end=30), bases="CAGGTGGATCATGAGGTCAGGAGTTCAAGA")
>>> p2: Primer = Primer(tm=37, penalty=0, span=Span(refname="chr1", start=61, end=93, strand=Strand.NEGATIVE), bases="CATGCCCAGCTAATTTTTTGTATTTTTAGTAGA")
>>> results_dict: dict[str, BwaResult] = detector.mappings_of(primers=[p1, p2])
>>> list(results_dict.keys())
['CAGGTGGATCATGAGGTCAGGAGTTCAAGA', 'CATGCCCAGCTAATTTTTTGTATTTTTAGTAGA']
>>> results = list(results_dict.values())
>>> results[0]
BwaResult(query=Query(id='CAGGTGGATCATGAGGTCAGGAGTTCAAGA', bases='CAGGTGGATCATGAGGTCAGGAGTTCAAGA'), hit_count=1, hits=[...])
>>> results[0].hits[0]
BwaHit(refname='chr1', start=1, negative=False, cigar=Cigar(elements=(CigarElement(length=30, operator=<CigarOp.M: (0, 'M', True, True)>),)), edits=0)
>>> results[1]
BwaResult(query=Query(id='CATGCCCAGCTAATTTTTTGTATTTTTAGTAGA', bases='CATGCCCAGCTAATTTTTTGTATTTTTAGTAGA'), hit_count=1, hits=[...])
>>> results[1].hits[0]
BwaHit(refname='chr1', start=61, negative=True, cigar=Cigar(elements=(CigarElement(length=33, operator=<CigarOp.M: (0, 'M', True, True)>),)), edits=0)

```

"""  # noqa: E501

import itertools
from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace
from pathlib import Path
from types import TracebackType
from typing import Optional
from typing import Self
from typing import TypeVar

from ordered_set import OrderedSet

from prymer.api.primer import Primer
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.offtarget.bwa import BwaAlnInteractive
from prymer.offtarget.bwa import BwaHit
from prymer.offtarget.bwa import BwaResult
from prymer.offtarget.bwa import Query

PrimerType = TypeVar("PrimerType", bound=Primer)


@dataclass(init=True, frozen=True)
class OffTargetResult:
    """Information obtained by running a single primer pair through the off-target detector.

    Attributes:
        primer_pair: the primer submitted
        passes: True if the primer pair passes all checks, False otherwise
        cached: True if this result is part of a cache, False otherwise.  This is useful for testing
        spans: the set of mappings of the primer pair to the genome or an empty list if mappings
            were not retained
        left_primer_spans: the list of mappings for the left primer, independent of the pair
            mappings, or an empty list
        right_primer_spans: the list of mappings for the right primer, independent of the pair
            mappings, or an empty list
    """

    primer_pair: PrimerPair
    passes: bool
    cached: bool = False
    spans: list[Span] = field(default_factory=list)
    left_primer_spans: list[Span] = field(default_factory=list)
    right_primer_spans: list[Span] = field(default_factory=list)


class OffTargetDetector:
    """A class for detecting off-target mappings of primers and primer pairs that uses a custom
    version of "bwa aln".

    The off-target detection is faster and more sensitive than traditional isPCR and in addition can
    correctly detect primers that are repetitive and contain many thousands or millions of mappings
    to the genome.

    Note that while this class invokes BWA with multiple threads, it is not itself thread-safe.
    Only one thread at a time should invoke methods on this class without external synchronization.
    """

    def __init__(
        self,
        ref: Path,
        max_primer_hits: int,
        max_primer_pair_hits: int,
        three_prime_region_length: int,
        max_mismatches_in_three_prime_region: int,
        max_mismatches: int,
        max_amplicon_size: int,
        cache_results: bool = True,
        threads: Optional[int] = None,
        keep_spans: bool = True,
        keep_primer_spans: bool = True,
        executable: str | Path = "bwa",
    ) -> None:
        """
        Args:
            ref: the reference genome fasta file (must be indexed with BWA)
            max_primer_hits: the maximum number of hits an individual primer can have in the genome
                before it is considered an invalid primer, and all primer pairs containing the
                primer failed.
            max_primer_pair_hits: the maximum number of amplicons a primer pair can make and be
                considered passing
            three_prime_region_length: the number of bases at the 3' end of the primer in which the
                parameter max_mismatches_in_three_prime_region is evaluated
            max_mismatches_in_three_prime_region: the maximum number of mismatches that are
                tolerated in the three prime region of each primer defined by
                three_prime_region_length
            max_mismatches: the maximum number of mismatches allowed in the full length primer
                (including any in the three prime region)
            max_amplicon_size: the maximum amplicon size to consider amplifiable
            cache_results: if True, cache results for faster re-querying
            threads: the number of threads to use when invoking bwa
            keep_spans: if True, [[OffTargetResult]] objects will have amplicon spans
                populated, otherwise not
            keep_primer_spans: if True, [[OffTargetResult]] objects will have left and right
                primer spans
            executable: string or Path representation of the `bwa` executable path
        """
        self._primer_cache: dict[str, BwaResult] = {}
        self._primer_pair_cache: dict[PrimerPair, OffTargetResult] = {}
        self._bwa = BwaAlnInteractive(
            executable=executable,
            ref=ref,
            reverse_complement=True,
            threads=threads,
            seed_length=three_prime_region_length,
            max_mismatches_in_seed=max_mismatches_in_three_prime_region,
            max_mismatches=max_mismatches,
            max_hits=max_primer_hits,
        )
        self._max_primer_hits = max_primer_hits
        self._max_primer_pair_hits: int = max_primer_pair_hits
        self._max_amplicon_size: int = max_amplicon_size
        self._cache_results: bool = cache_results
        self._keep_spans: bool = keep_spans
        self._keep_primer_spans: bool = keep_primer_spans

    def filter(self, primers: list[PrimerType]) -> list[PrimerType]:
        """Filters an iterable of Primers to return only those that have less than
        `max_primer_hits` mappings to the genome."""
        results: dict[str, BwaResult] = self.mappings_of(primers)
        return [
            primer for primer in primers if results[primer.bases].hit_count <= self._max_primer_hits
        ]

    def check_one(self, primer_pair: PrimerPair) -> OffTargetResult:
        """Checks a PrimerPair for off-target sites in the genome at which it might amplify."""
        result: dict[PrimerPair, OffTargetResult] = self.check_all([primer_pair])
        return result[primer_pair]

    def check_all(self, primer_pairs: list[PrimerPair]) -> dict[PrimerPair, OffTargetResult]:
        """Checks a collection of primer pairs for off-target sites, returning a dictionary of
        `PrimerPair`s to `OffTargetResult`.

        Returns:
            a Map containing all given primer pairs as keys with the values being the result of
            off-target checking.
        """

        primer_pair_results: dict[PrimerPair, OffTargetResult] = {}
        result: OffTargetResult

        # Get the primer pairs to map.  If the primer pair is found in the cache, use that
        primer_pairs_to_map: list[PrimerPair] = []
        if not self._cache_results:
            primer_pairs_to_map = primer_pairs
        else:
            for primer_pair in primer_pairs:
                match self._primer_pair_cache.get(primer_pair, None):
                    case None:
                        primer_pairs_to_map.append(primer_pair)  # Map it!
                    case result:
                        primer_pair_results[primer_pair] = result

        # If there are no primer pairs to map, return the results
        if len(primer_pairs_to_map) == 0:
            return primer_pair_results

        # Get mappings of all the primers
        primers = [primer for primer_pair in primer_pairs for primer in primer_pair]
        hits_by_primer = self.mappings_of(primers)

        for primer_pair in primer_pairs:
            primer_pair_results[primer_pair] = self._build_off_target_result(
                primer_pair=primer_pair, hits_by_primer=hits_by_primer
            )

        return primer_pair_results

    def _build_off_target_result(
        self, primer_pair: PrimerPair, hits_by_primer: dict[str, BwaResult]
    ) -> OffTargetResult:
        """Builds an `OffTargetResult` for the given `PrimerPair`.

        If there are too many primer hits for either the left or right primer, set `passes` to
        `False` on the returned `OffTargetResult`, otherwise generate all possible amplicons with
        size less than the given maximum to generate off target results.

        Args:
            primer_pair: the primer pair from which the off-target result is built.
            hits_by_primer: mapping of primer bases to `BwaHits`.
        """
        result: OffTargetResult

        # Get the mappings for the left primer and right primer respectively
        p1: BwaResult = hits_by_primer[primer_pair.left_primer.bases]
        p2: BwaResult = hits_by_primer[primer_pair.right_primer.bases]

        # Get all possible amplicons from the left_primer_mappings and right_primer_mappings
        # primer hits, filtering if there are too many for either
        if p1.hit_count > self._max_primer_hits or p2.hit_count > self._max_primer_hits:
            result = OffTargetResult(primer_pair=primer_pair, passes=False)
        else:
            amplicons = self._to_amplicons(p1.hits, p2.hits, self._max_amplicon_size)
            result = OffTargetResult(
                primer_pair=primer_pair,
                passes=0 < len(amplicons) <= self._max_primer_pair_hits,
                spans=amplicons if self._keep_spans else [],
                left_primer_spans=(
                    [self._hit_to_span(h) for h in p1.hits] if self._keep_primer_spans else []
                ),
                right_primer_spans=(
                    [self._hit_to_span(h) for h in p2.hits] if self._keep_primer_spans else []
                ),
            )

        if self._cache_results:
            self._primer_pair_cache[primer_pair] = replace(result, cached=True)

        return result

    def mappings_of(self, primers: list[PrimerType]) -> dict[str, BwaResult]:
        """Function to take a set of primers and map any that are not cached, and return mappings
        for all of them.  Note: the genomics sequence of the returned mappings are on the opposite
        strand of that of the strand of the primer.  I.e. we map the complementary bases (reversed)
        to that of the primer."""

        primers_to_map: list[PrimerType]
        if not self._cache_results:
            primers_to_map = primers
        else:
            primers_to_map = [
                primer for primer in list(OrderedSet(primers)) if primer not in self._primer_cache
            ]

        # Build the unique list of queries to map with BWA
        queries: list[Query] = [
            Query(id=primer.bases, bases=primer.bases) for primer in primers_to_map
        ]

        # Map the queries with BWA
        hits_by_primer: dict[str, BwaResult] = {
            result.query.id: result for result in self._bwa.map_all(queries)
        }

        # Cache the results, if desired, and get the hits by primer _for all_ primers. If not
        # caching, then hits_by_primer already contains all the primers.
        if self._cache_results:
            self._primer_cache.update(hits_by_primer)
            hits_by_primer = {primer.bases: self._primer_cache[primer.bases] for primer in primers}

        return hits_by_primer

    @staticmethod
    def _to_amplicons(lefts: list[BwaHit], rights: list[BwaHit], max_len: int) -> list[Span]:
        """Takes a set of hits for one or more left primers and right primers and constructs
        amplicon mappings anywhere a left primer hit and a right primer hit align in F/R
        orientation up to `maxLen` apart on the same reference.  Primers may not overlap.
        """
        amplicons: list[Span] = []
        for h1, h2 in itertools.product(lefts, rights):
            if h1.negative == h2.negative or h1.refname != h2.refname:  # not F/R orientation
                continue

            plus, minus = (h2, h1) if h1.negative else (h1, h2)
            if minus.start > plus.end and (minus.end - plus.start + 1) <= max_len:
                amplicons.append(Span(refname=plus.refname, start=plus.start, end=minus.end))

        return amplicons

    @staticmethod
    def _hit_to_span(hit: BwaHit) -> Span:
        """Converts a Bwa Hit object to a Span."""
        return Span(refname=hit.refname, start=hit.start, end=hit.end)

    def close(self) -> None:
        self._bwa.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> None:
        """Gracefully terminates any running subprocesses."""
        self.close()
