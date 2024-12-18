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
>>> left_primer =  Oligo(bases="AAAAA", tm=37, penalty=0, span=Span("chr1", start=67, end=71))
>>> right_primer = Oligo(bases="TTTTT", tm=37, penalty=0, span=Span("chr1", start=75, end=79, strand=Strand.NEGATIVE))
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
>>> p1: Oligo = Oligo(tm=37, penalty=0, span=Span(refname="chr1", start=1, end=30), bases="CAGGTGGATCATGAGGTCAGGAGTTCAAGA")
>>> p2: Oligo = Oligo(tm=37, penalty=0, span=Span(refname="chr1", start=61, end=93, strand=Strand.NEGATIVE), bases="CATGCCCAGCTAATTTTTTGTATTTTTAGTAGA")
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
from contextlib import AbstractContextManager
from dataclasses import dataclass
from dataclasses import field
from dataclasses import replace
from pathlib import Path
from types import TracebackType
from typing import Optional
from typing import Self
from typing import TypeVar

from ordered_set import OrderedSet

from prymer.api.oligo import Oligo
from prymer.api.primer_pair import PrimerPair
from prymer.api.span import Span
from prymer.offtarget.bwa import BWA_EXECUTABLE_NAME
from prymer.offtarget.bwa import BwaAlnInteractive
from prymer.offtarget.bwa import BwaHit
from prymer.offtarget.bwa import BwaResult
from prymer.offtarget.bwa import Query

PrimerType = TypeVar("PrimerType", bound=Oligo)


@dataclass(init=True, frozen=True)
class OffTargetResult:
    """Information obtained by running a single primer pair through the off-target detector.

    Attributes:
        primer_pair: the primer pair submitted
        passes: True if neither primer exceeds the specified maximum number of hits, and the primer
            pair would generate an acceptable number of amplicons in the reference genome (i.e.
            `min_primer_pair_hits <= num_amplicons <= max_primer_pair_hits`). False otherwise.
        cached: True if this result is part of a cache, False otherwise.  This is useful for testing
        spans: The set of mappings of the primer pair to the genome (i.e., the region spanned by the
            inferred amplicon). This list will be empty if the generating [[OffTargetDetector]] was
            constructed with `keep_spans=False`.
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


class OffTargetDetector(AbstractContextManager):
    """A class for detecting off-target mappings of primers and primer pairs that uses a custom
    version of "bwa aln" named "bwa-aln-interactive".

    `OffTargetDetector` uses a [custom, interactive
    version](https://github.com/fulcrumgenomics/bwa-aln-interactive/) of `bwa aln` to perform
    off-target detection. This approach is faster and more sensitive than traditional isPCR and in
    addition can correctly detect primers that are repetitive and contain many thousands or millions
    of mappings to the genome.

    Note that while this class invokes BWA with multiple threads, it is not itself thread-safe.
    Only one thread at a time should invoke methods on this class without external synchronization.

    Off-target detection may be performed for individual primers and/or primer pairs.

    To detect off-target hits for individual primers, use the `OffTargetDetector.filters()` method,
    which will remove any primers that have more than the specified number of maximum hits against
    the reference.

    To detect off-target amplification of primer pairs, use the `OffTargetDetector.check_one()` or
    `OffTargetDetector.check_all()` methods. These methods screen the individual primers in each
    pair for off-target hits, and verify that the possible amplicons inferred by the primers'
    alignments does not exceed the specified maximum number of primer pair hits.
    """

    def __init__(
        self,
        ref: Path,
        max_primer_hits: int,
        max_primer_pair_hits: int,
        three_prime_region_length: int,
        max_mismatches_in_three_prime_region: int,
        max_mismatches: int,
        max_gap_opens: int = 0,
        max_gap_extends: int = -1,
        max_amplicon_size: int = 1000,
        min_primer_pair_hits: int = 1,
        cache_results: bool = True,
        threads: Optional[int] = None,
        keep_spans: bool = True,
        keep_primer_spans: bool = True,
        executable: str | Path = BWA_EXECUTABLE_NAME,
    ) -> None:
        """
        Initialize an [[OffTargetDetector]].

        This class accepts a large number of parameters to control the behavior of the off-target
        detection. The parameters are used in the various aspects of off-target detection as
        follows:

        1. Alignment (off-target hit detection): `ref`, `executable`, `threads`,
           `three_prime_region_length`, `max_mismatches_in_three_prime_region`, `max_mismatches`,
           and `max_primer_hits`.
        2. Filtering of individual primers: `max_primer_hits`.
        3. Checking of primer pairs: `max_primer_hits`, `min_primer_pair_hits`,
           `max_primer_pair_hits`, and `max_amplicon_size`.

        Args:
            ref: the reference genome fasta file (must be indexed with BWA)
            max_primer_hits: the maximum number of hits an individual primer can have in the genome
                before it is considered an invalid primer, and all primer pairs containing the
                primer failed.
            max_primer_pair_hits: the maximum number of amplicons a primer pair can make and be
                considered passing
            min_primer_pair_hits: The minimum number of amplicons a primer pair can make and be
                considered passing. (In most cases, this is the number of amplicons a primer pair is
                expected to generate.) The default is 1, which is appropriate when the primer pair
                is being evaluated for off-target hits against the same reference genome from which
                the primers were generated. If the primer pair was generated from a different
                reference sequence, it may be appropriate to set this value to 0.
            three_prime_region_length: the number of bases at the 3' end of the primer in which the
                parameter max_mismatches_in_three_prime_region is evaluated
            max_mismatches_in_three_prime_region: the maximum number of mismatches that are
                tolerated in the three prime region of each primer defined by
                three_prime_region_length
            max_mismatches: the maximum number of mismatches allowed in the full length primer
                (including any in the three prime region)
            max_gap_opens: the maximum number of gaps (insertions or deletions) allowable in an
                alignment of a oligo to the reference
            max_gap_extends: the maximum number of gap extensions allowed; extending a gap
                beyond a single base costs 1 gap extension.  Can be set to -1 to allow
                unlimited extensions up to max diffs (aka max mismatches), while disallowing
                "long gaps".
            max_amplicon_size: the maximum amplicon size to consider amplifiable
            cache_results: if True, cache results for faster re-querying
            threads: the number of threads to use when invoking bwa
            keep_spans: if True, [[OffTargetResult]] objects will be reported with amplicon spans
                populated, otherwise not
            keep_primer_spans: if True, [[OffTargetResult]] objects will be reported with left and
                right primer spans
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
            max_gap_opens=max_gap_opens,
            max_gap_extensions=max_gap_extends,
            max_hits=max_primer_hits,
        )
        self._max_primer_hits = max_primer_hits
        self._max_primer_pair_hits: int = max_primer_pair_hits
        self._min_primer_pair_hits: int = min_primer_pair_hits
        self._max_amplicon_size: int = max_amplicon_size
        self._cache_results: bool = cache_results
        self._keep_spans: bool = keep_spans
        self._keep_primer_spans: bool = keep_primer_spans

    def filter(self, primers: list[PrimerType]) -> list[PrimerType]:
        """
        Remove primers that have more than `max_primer_hits` mappings to the genome.

        This method maps each primer to the specified reference with `bwa aln` to search for
        off-target hits. Note that when the reference includes the sequence from which the primers
        were designed, the on-target hit will be included in the hit count. `max_primer_hits` should
        be set to at least 1 in this case.

        Args:
            primers: A list of primers to filter.

        Returns:
            The input primers, filtered to remove any primers with more than `max_primer_hits`
            mappings to the genome.
        """
        results: dict[str, BwaResult] = self.mappings_of(primers)
        return [
            primer for primer in primers if results[primer.bases].hit_count <= self._max_primer_hits
        ]

    def check_one(self, primer_pair: PrimerPair) -> OffTargetResult:
        """
        Check an individual primer pair for off-target amplification.

        See `OffTargetDetector.check_all()` for details.

        Args:
            primer_pair: The primer pair to check.

        Returns:
            The results of the off-target detection.
            See `OffTargetResult` for details regarding the attributes included in the result.
        """
        result: dict[PrimerPair, OffTargetResult] = self.check_all([primer_pair])
        return result[primer_pair]

    def check_all(self, primer_pairs: list[PrimerPair]) -> dict[PrimerPair, OffTargetResult]:
        """
        Check a collection of primer pairs for off-target amplification.

        This method maps each primer to the specified reference with `bwa aln` to search for
        off-target hits. Possible amplicons are identified from the pairwise combinations of these
        alignments, up to the specified `max_amplicon_size`.

        Primer pairs are marked as passing if both of the following are true:
        1. Each primer has no more than `max_primer_hits` alignments to the genome.
        2. The size of the set of possible amplicons does not exceed `max_primer_pair_hits`, and is
           at least `min_primer_pair_hits`.

        Args:
            primer_pairs: The primer pairs to check.

        Returns:
            A dictionary mapping each checked primer pair to its off-target detection results.
            See `OffTargetResult` for details regarding the attributes included in each result.
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
                passes=self._min_primer_pair_hits <= len(amplicons) <= self._max_primer_pair_hits,
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
        """
        Map primers to the reference genome using `bwa aln`.

        Alignment results may be optionally cached for efficiency. Set `cache_results` to `True`
        when instantiating the [[OffTargetDetector]] to enable caching.

        Any primers without cached results, or all primers when `cache_results` is `False`, will be
        mapped to the reference genome using `bwa aln` and the specified alignment parameters.

        **Note**: The reverse complement of each primer sequence is used for mapping. The `query`
        sequence in each `BwaResult` will match the input primer sequence, as will the sequences
        used as keys in the output dictionary. However, the coordinates reported in each `BwaHit`
        associated with a result will correspond to complementary sequence.

        Args:
            primers: A list of primers to map.

        Returns:
            A dictionary mapping each primer's sequence to its alignment results.
            See `BwaResult` for details regarding the attributes included in each result.
        """

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
        super().__exit__(exc_type, exc_value, traceback)
        self.close()
