"""
# Methods for merging intervals into "clusters"

This module contains utility functions related to clustering intervals into larger intervals.

There is currently one public method available:

- [`cluster_intervals()`][prymer.api.clustering.cluster_intervals] -- clusters a list of
  intervals into a set of intervals that cover the input intervals, and are not larger than an input
  `max_size`. each original interval will be wholly covered by at least one cluster.  The `name`
  attribute of each cluster will be set, and the original intervals are returned with a `name`
  attribute that matches that of a cluster that wholly contains it.

## Examples

```python
>>> from prymer.api.clustering import cluster_intervals
>>> from pybedlite.overlap_detector import Interval
>>> intervals = [Interval("chr1", 1, 2), Interval("chr1", 3, 4)]
>>> cluster_intervals(intervals, 10)
ClusteredIntervals(clusters=[Interval(refname='chr1', start=1, end=4, negative=False, name='chr1:1-4')], intervals=[Interval(refname='chr1', start=1, end=2, negative=False, name='chr1:1-4'), Interval(refname='chr1', start=3, end=4, negative=False, name='chr1:1-4')])
>>> cluster_intervals(intervals, 2)
ClusteredIntervals(clusters=[Interval(refname='chr1', start=1, end=2, negative=False, name='chr1:1-2'), Interval(refname='chr1', start=3, end=4, negative=False, name='chr1:3-4')], intervals=[Interval(refname='chr1', start=1, end=2, negative=False, name='chr1:1-2'), Interval(refname='chr1', start=3, end=4, negative=False, name='chr1:3-4')])

```
"""  # noqa: E501

import itertools
from dataclasses import dataclass
from typing import Dict

from attr import evolve
from pybedlite.overlap_detector import Interval
from pybedlite.overlap_detector import OverlapDetector

from prymer.api.coordmath import get_locus_string
from prymer.api.coordmath import require_same_refname


@dataclass(frozen=True, init=True, slots=True)
class ClusteredIntervals:
    """
    The list of clusters (intervals) and the original source intervals.  The source intervals must
    have the name corresponding to the cluster to which the source interval belongs. Each cluster
    must envelop ("wholly contain") the intervals associated with the cluster.

    Attributes:
        clusters: the clusters that wholly contain one or more source intervals.
        intervals: the source intervals, with name corresponding to the name of the associated
            cluster.
    """

    clusters: list[Interval]
    intervals: list[Interval]

    def __post_init__(self) -> None:
        cluster_names: set[str] = {c.name for c in self.clusters}
        interval_names: set[str] = {i.name for i in self.intervals}
        # Check that the names of clusters are unique
        if len(self.clusters) != len(cluster_names):
            raise ValueError("Cluster names are not unique")
        # Check that every interval has the name of one cluster
        for interval in self.intervals:
            if interval.name not in cluster_names:
                raise ValueError(f"Interval does not belong to a cluster: {interval}")
        # Check that every cluster has at least one interval associated with the cluster
        if len(interval_names) != len(cluster_names):
            raise ValueError("Cluster and interval names differ.")


def _sort_key(interval: Interval) -> tuple[int, int]:
    """Returns the sort key to use when sorting intervals from the same reference

    Note that this method assumes that the intervals have the same reference name, but this is NOT
    checked.  Using it for intervals from different reference will result in undefined behavior.

    Args:
        interval: The interval to get the sort key for.
    Returns:
        A tuple (start, end) of the interval.
    """
    return interval.start, interval.end


def _cluster_in_contig(intervals: list[Interval], max_size: int) -> ClusteredIntervals:
    """
    Cluster a list of intervals (all from one reference) into intervals that overlap the given
    intervals and are not larger than `max_size`.

    Implements a greedy algorithm for hierarchical clustering, merging subsequent intervals
    (from a sorted list) as long as the maximal size is respected.
    Each "cluster" is replaced by an interval that spans it, and the algorithm terminates
    when it can no longer merge anything without creating a cluster that is larger than `max_size`.

    Args:
        intervals: The intervals to cluster.
        max_size: The maximum size (in bp) of the resulting clusters.

    Returns:
        A named tuple (`clusters`, `intervals`), where `clusters` contains one (named) interval per
        cluster, defining the region spanned by the cluster, and `intervals` contains the original
        set of intervals, each adorned with a `name` that agrees with that of a cluster in
        `clusters` that wholly contains it.

    Raises:
        ValueError: If any of the input intervals are larger than `max_size`.
        ValueError: If the input intervals contain between them more than one refname.
    """
    if len(intervals) == 0:
        return ClusteredIntervals(clusters=[], intervals=[])

    require_same_refname(*intervals)
    max_found = max(i.length() for i in intervals)
    if max_found > max_size:
        raise ValueError(
            f"Intervals provided must be less than {max_size}, " f"but found size {max_found}"
        )

    sorted_intervals = sorted(intervals, key=_sort_key)
    # at each step, check if can "legally" merge the current cluster with the next interval,
    # and if so, update the current cluster, if not, add the current cluster to the list and
    # start a new cluster with the current interval.
    clusters: list[Interval] = []

    curr_cluster = sorted_intervals[0]
    for interval in sorted_intervals[1:]:
        new_cluster = _convex_hull(curr_cluster, interval)

        if new_cluster.length() <= max_size:
            curr_cluster = new_cluster
        else:
            clusters.append(evolve(curr_cluster, name=get_locus_string(curr_cluster)))
            curr_cluster = interval

    clusters.append(evolve(curr_cluster, name=get_locus_string(curr_cluster)))

    # for each original interval find one cluster that it is contained in.
    detector: OverlapDetector[Interval] = OverlapDetector()
    detector.add_all(clusters)

    enclosing_clusters = [detector.get_enclosing_intervals(si).pop() for si in sorted_intervals]
    ann_intervals: list[Interval] = [
        evolve(interval, name=enc_cluster.name)
        for interval, enc_cluster in zip(sorted_intervals, enclosing_clusters, strict=True)
    ]

    return ClusteredIntervals(clusters=clusters, intervals=ann_intervals)


def cluster_intervals(
    intervals: list[Interval],
    max_size: int,
) -> ClusteredIntervals:
    """
    Cluster a list of intervals into intervals that overlap the given
    intervals and are not larger than `max_size`.

    Implements a greedy algorithm for hierarchical clustering, merging subsequent intervals
    (from a sorted list) as long as the maximal size is respected.
    Each "cluster" is replaced by an interval that spans it, and the algorithm terminates
    when it can no longer merge anything without creating a cluster that is larger than `max_size`.

    Args:
        intervals: The intervals to cluster.
        max_size: The maximum size (in bp) of the resulting clusters.

    Returns:
        A named tuple (`clusters`, `intervals`), where `clusters` contains one (named) interval per
        cluster, defining the region spanned by the cluster, and `intervals` contains the original
        set of intervals, each adorned with a `name` that agrees with that of a cluster in
        `clusters` that wholly contains it.

    Raises:
        ValueError: If any of the input intervals are larger than `max_size`.
    """

    intervals_by_refname: Dict[str, list[Interval]] = {
        key: list(group) for key, group in itertools.groupby(intervals, key=lambda x: x.refname)
    }
    # import ipdb; ipdb.set_trace()
    # iterate over the refnames and call _cluster_in_contig() per refname.

    # need to store this in a list since we iterate over it twice.
    per_refname_clusters_and_intervals = [
        x
        for x in (
            _cluster_in_contig(intervals_in_refname, max_size)
            for intervals_in_refname in intervals_by_refname.values()
        )
    ]

    # flatten the clusters and intervals
    clusters = [
        cluster
        for refname_result in per_refname_clusters_and_intervals
        for cluster in refname_result.clusters
    ]

    intervals_out = [
        interval
        for refname_results in per_refname_clusters_and_intervals
        for interval in refname_results.intervals
    ]

    return ClusteredIntervals(clusters=clusters, intervals=intervals_out)


def _convex_hull(interval_a: Interval, interval_b: Interval) -> Interval:
    """
    Get the convex hull of two intervals.

    Args:
        interval_a: The first interval.
        interval_b: The second interval.
    Returns:
        The convex hull of the two intervals (a new interval).

    Raises:
        ValueError: If the intervals do not have the same refname.
    """
    require_same_refname(interval_a, interval_b)
    return evolve(
        interval_a,
        start=min(interval_a.start, interval_b.start),
        end=max(interval_a.end, interval_b.end),
        name="",
    )
