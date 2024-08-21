import io

import pytest
from pybedlite.bed_source import BedSource
from pybedlite.overlap_detector import Interval

from prymer.api.clustering import _cluster_in_contig
from prymer.api.clustering import cluster_intervals

with BedSource(
    io.StringIO("""
chr1\t3\t6\tone\t0\t+
chr1\t8\t9\ttwo\t0\t+
chr1\t5\t7\tthree\t0\t+
""")
) as bed:
    d = [Interval.from_bedrecord(i) for i in bed]

with BedSource(
    io.StringIO("""
chr1\t3\t6\tone\t0\t+
chr1\t8\t9\ttwo\t0\t+
chr1\t5\t7\tthree\t0\t+
chr2\t3\t6\tone\t0\t+
chr2\t8\t9\ttwo\t0\t+
chr2\t5\t7\tthree\t0\t+
""")
) as bed:
    both = [Interval.from_bedrecord(i) for i in bed]


def test_clustering_when_no_overlaps() -> None:
    no_overlap = [d[i] for i in [0, 1]]
    result = cluster_intervals(no_overlap, 10)

    # print(clusters)

    assert len(result.clusters) == 1
    assert len(result.intervals) == 2
    assert len({i.name for i in result.clusters}) == 1


def test_cluster_in_contig() -> None:
    result = _cluster_in_contig(d, 10)

    assert len(result.clusters) == 1
    assert len(result.intervals) == len(d)


def test_cluster_empty() -> None:
    result = _cluster_in_contig([], 10)

    assert len(result.clusters) == 0
    assert len(result.intervals) == 0


def test_cluster_in_contig_small_distance() -> None:
    result = _cluster_in_contig(d, 3)
    assert len(result.clusters) == 3
    assert len(result.intervals) == len(d)
    assert len({i.name for i in result.intervals}) == 3


def test_cluster_fails_on_tiny_distance() -> None:
    with pytest.raises(ValueError):
        _cluster_in_contig(d, 2)


def test_cluster_multi_contig() -> None:
    result = cluster_intervals(both, 10)

    assert len(result.clusters) == 2
    assert len({i.refname for i in result.clusters}) == 2
    assert len(result.intervals) == 2 * len(d)
    assert len({i.name for i in result.intervals}) == 2


def test_cluster_multi_contig_small_distance() -> None:
    result = cluster_intervals(both, 3)

    assert max(i.length() for i in result.clusters) <= 3
    assert len(result.clusters) == 2 * 3
    assert len({i.refname for i in result.clusters}) == 2
    assert len(result.intervals) == 2 * len(d)
    assert len({i.name for i in result.intervals}) == 3 * 2


def test_cluster_multi_contig_fails_on_tiny_distance() -> None:
    with pytest.raises(ValueError):
        cluster_intervals(both, 2)


def test_max_size_respected() -> None:
    n = 100
    max_size = 30
    interval_size = 5
    intervals = [Interval("chr1", i, i + interval_size) for i in range(n)]

    result = cluster_intervals(intervals, max_size)

    assert max(i.length() for i in result.clusters) <= max_size
