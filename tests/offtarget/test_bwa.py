import shutil
from dataclasses import replace
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import TypeAlias
from typing import cast

import pytest
from fgpyo.sam import Cigar
from fgpyo.sam.builder import SamBuilder

from prymer.offtarget.bwa import BWA_AUX_EXTENSIONS
from prymer.offtarget.bwa import BwaAlnInteractive
from prymer.offtarget.bwa import BwaHit
from prymer.offtarget.bwa import Query

SamHeaderType: TypeAlias = dict[str, dict[str, str] | list[dict[str, str]]]
"""A type alias for a SAM sequence header dictionary."""


@pytest.mark.parametrize("bases", [None, ""])
def test_query_no_bases(bases: str | None) -> None:
    with pytest.raises(ValueError, match="No bases were provided to query"):
        Query(id="foo", bases=bases)


def test_bwa_hit_mismatches_indels_gt_edits() -> None:
    hit = BwaHit.build("chr1", 1081, False, "20M5I30M", 0, revcomp=False)
    with pytest.raises(ValueError, match="indel_sum"):
        assert 0 == hit.mismatches


def test_bwa_hit_mismatches() -> None:
    hit = BwaHit.build("chr1", 1081, False, "20M5I30M", 5, revcomp=False)
    assert 0 == hit.mismatches
    hit = BwaHit.build("chr1", 1081, False, "22M3I30M", 5, revcomp=False)
    assert 2 == hit.mismatches
    hit = BwaHit.build("chr1", 1081, False, "55M", 5, revcomp=False)
    assert 5 == hit.mismatches


def test_missing_index_files(genome_ref: Path) -> None:
    # copy to temp dir
    with NamedTemporaryFile(suffix=".fa", mode="w", delete=True) as temp_file:
        ref_fasta = Path(temp_file.name)
        shutil.copy(genome_ref, ref_fasta)

        # no files exist
        with pytest.raises(FileNotFoundError, match="BWA index files do not exist"):
            BwaAlnInteractive(ref=ref_fasta, max_hits=1)

        # all but one missing
        for aux_ext in BWA_AUX_EXTENSIONS[1:]:
            aux_path = Path(f"{ref_fasta}{aux_ext}")
            with aux_path.open("w"):
                pass
        with pytest.raises(FileNotFoundError, match="BWA index file does not exist"):
            BwaAlnInteractive(ref=ref_fasta, max_hits=1)


def test_hit_build_rc() -> None:
    hit_f = BwaHit.build("chr1", 1081, False, "20M5I30M", 5, revcomp=False)
    hit_r = BwaHit.build("chr1", 1081, False, "20M5I30M", 5, revcomp=True)
    # things that stay the same
    assert hit_f.refname == hit_r.refname
    assert hit_f.start == hit_r.start
    assert hit_f.end == hit_r.end
    assert hit_f.mismatches == hit_r.mismatches
    assert hit_f.edits == hit_r.edits
    assert hit_f.cigar == hit_r.cigar
    # things that change
    assert hit_f.negative != hit_r.negative
    assert hit_r.negative is True


def test_header_is_properly_constructed(ref_fasta: Path) -> None:
    """Tests that bwa will return a properly constructed header."""
    with BwaAlnInteractive(ref=ref_fasta, max_hits=1) as bwa:
        header: SamHeaderType = bwa.header.to_dict()
        assert set(header.keys()) == {"HD", "SQ", "PG"}
        assert header["HD"] == {"GO": "query", "SO": "unsorted", "VN": "1.5"}
        assert header["SQ"] == [{"LN": 10001, "SN": "chr1"}]
        assert len(header["PG"]) == 1
        program_group: dict[str, str] = cast(list[dict[str, str]], header["PG"])[0]
        program_group["ID"] = "bwa"
        program_group["PN"] = "bwa"


def test_map_one_uniquely_mapped(ref_fasta: Path) -> None:
    """Tests that bwa maps one hit when a query uniquely maps."""
    query = Query(bases="TCTACTAAAAATACAAAAAATTAGCTGGGCATGATGGCATGCACCTGTAATCCCGCTACT", id="NA")
    with BwaAlnInteractive(ref=ref_fasta, max_hits=1) as bwa:
        result = bwa.map_one(query=query.bases, id=query.id)
        assert result.hit_count == 1
        assert result.hits[0].refname == "chr1"
        assert result.hits[0].start == 61
        assert result.hits[0].negative is False
        assert f"{result.hits[0].cigar}" == "60M"
        assert result.query == query


def test_map_one_unmapped(ref_fasta: Path) -> None:
    """Tests that bwa returns an unmapped alignment.  The hit count should be zero and the list
    of hits empty."""
    with BwaAlnInteractive(ref=ref_fasta, max_hits=1) as bwa:
        query = Query(bases="A" * 50, id="NA")
        result = bwa.map_one(query=query.bases, id=query.id)
        assert result.hit_count == 0
        assert len(result.hits) == 0
        assert result.query == query


def test_map_one_multi_mapped_max_hits_one(ref_fasta: Path) -> None:
    """Tests that a query that returns too many hits (>max_hits) returns the number of hits but
    not the list of hits themselves."""
    with BwaAlnInteractive(ref=ref_fasta, max_hits=1) as bwa:
        query = Query(bases="A" * 5, id="NA")
        result = bwa.map_one(query=query.bases, id=query.id)
        assert result.hit_count == 7508
        assert len(result.hits) == 0  # hit_count > max_hits
        assert result.query == query


def test_map_one_multi_mapped_max_hits_many(ref_fasta: Path) -> None:
    """Tests a query that aligns to many locations, but fewer than max_hits, returns the number of
    hits and the hits themselves"""
    with BwaAlnInteractive(ref=ref_fasta, max_hits=10000) as bwa:
        query = Query(bases="A" * 5, id="NA")
        result = bwa.map_one(query=query.bases, id=query.id)
        assert result.hit_count == 7504
        assert len(result.hits) == 7504  # hit_count <= max_hits
        assert result.query == query


def test_map_all(ref_fasta: Path) -> None:
    """Tests aligning multiple queries."""
    with BwaAlnInteractive(ref=ref_fasta, max_hits=10000) as bwa:
        # empty queries
        assert bwa.map_all([]) == []

        # many queries
        queries = []
        for length in range(5, 101):
            queries.append(Query(bases="A" * length, id="NA"))
        results = bwa.map_all(queries)
        assert len(results) == len(queries)
        assert all(results[i].query == queries[i] for i in range(len(results)))
        assert results[0].hit_count == 7504
        assert len(results[0].hits) == 7504  # hit_count <= max_hits
        assert results[0].query == queries[0]
        assert results[-1].hit_count == 0
        assert len(results[-1].hits) == 0
        assert results[-1].query == queries[-1]


_PERFECT_BASES: str = "AGTGATGCTAAGGGTCAAATAAGTCACCAGCAAATACACAGCACACATCTCATGATGTGC"
"""Bases for perfect matching hit to the miniref.fa"""

_PERFECT_HIT: BwaHit = BwaHit.build("chr1", 1081, False, "60M", 0)
"""Perfect matching hit to the miniref.fa"""


# TODO: # of indels
@pytest.mark.parametrize(
    "mismatches, end, bases, hit",
    [
        # Perfect-matching hit
        (0, 1140, _PERFECT_BASES, _PERFECT_HIT),
        # One mismatch hit
        (1, 1140, _PERFECT_BASES[:30] + "N" + _PERFECT_BASES[31:], replace(_PERFECT_HIT, edits=1)),
        # Five mismatch hit
        (
            5,
            1140,
            _PERFECT_BASES[:30] + "NNNNN" + _PERFECT_BASES[35:],
            replace(_PERFECT_HIT, edits=5),
        ),
        # Deletion
        (
            0,
            1140,
            _PERFECT_BASES[:31] + _PERFECT_BASES[36:],
            replace(_PERFECT_HIT, edits=5, cigar=Cigar.from_cigarstring("31M5D24M")),
        ),
        # Insertion
        (
            0,
            1140,
            _PERFECT_BASES[:31] + "TTTTT" + _PERFECT_BASES[31:],
            replace(_PERFECT_HIT, edits=5, cigar=Cigar.from_cigarstring("31M5I29M")),
        ),
    ],
)
@pytest.mark.parametrize("reverse_complement", [True, False])
def test_map_single_hit(
    ref_fasta: Path, mismatches: int, end: int, bases: str, hit: BwaHit, reverse_complement: bool
) -> None:
    """Tests that bwa maps one hit when a query uniquely maps.  Checks for the hit's edit
    and mismatches properties."""
    with BwaAlnInteractive(
        ref=ref_fasta,
        max_hits=1,
        max_mismatches=5,
        max_gap_opens=1,
        reverse_complement=reverse_complement,
    ) as bwa:
        query = Query(bases=bases, id=f"{hit}")
        result = bwa.map_one(query=query.bases, id=query.id)
        assert result.hit_count == 1, "hit_count"
        assert result.query == query, "query"
        assert result.hits[0].refname == hit.refname, "chr"
        assert result.hits[0].start == hit.start, "start"
        assert result.hits[0].negative == hit.negative, "negative"
        assert result.hits[0].cigar == hit.cigar, "cigar"
        assert result.hits[0].edits == hit.edits, "edits"
        assert result.hits[0].end == end, "end"
        assert result.hits[0].mismatches == mismatches, "mismatches"
        assert result.hits[0].mismatches == hit.mismatches, "hit.mismatches"


@pytest.mark.parametrize("reverse_complement", [True, False])
def test_map_multi_hit(ref_fasta: Path, reverse_complement: bool) -> None:
    mismatches: int = 0
    expected_hits: list[BwaHit] = [
        BwaHit.build("chr1", 2090, False, "31M", 0),
        BwaHit.build("chr1", 5825, False, "31M", 0),
        BwaHit.build("chr1", 8701, False, "31M", 0),
    ]
    bases: str = "AAGTGCTGGGATTACAGGCATGAGCCACCAC"
    with BwaAlnInteractive(
        ref=ref_fasta,
        max_hits=len(expected_hits),
        max_mismatches=mismatches,
        max_gap_opens=0,
        reverse_complement=reverse_complement,
    ) as bwa:
        query = Query(bases=bases, id="test")
        actual = bwa.map_one(query=query.bases, id=query.id)
        assert actual.hit_count == 3, "hit_count"
        assert actual.query == query, "query"
        actual_hits = sorted(actual.hits, key=lambda hit: (hit.refname, hit.start))
        assert len(actual_hits) == len(expected_hits)
        for actual_hit, expected_hit in zip(actual_hits, expected_hits, strict=True):
            assert actual_hit.refname == expected_hit.refname, "chr"
            assert actual_hit.start == expected_hit.start, "start"
            assert actual_hit.negative == expected_hit.negative, "negative"
            assert actual_hit.cigar == expected_hit.cigar, "cigar"
            assert actual_hit.edits == expected_hit.edits, "edits"
            assert actual_hit.end == expected_hit.start + 30, "end"
            assert actual_hit.mismatches == mismatches, "mismatches"
            assert actual_hit.mismatches == expected_hit.mismatches, "hit.mismatches"


def test_to_result_out_of_order(ref_fasta: Path) -> None:
    with BwaAlnInteractive(ref=ref_fasta, max_hits=1) as bwa:
        query = Query(bases="GATTACA", id="foo")
        rec = SamBuilder().add_single(name="bar")
        with pytest.raises(ValueError, match="Query and Results are out of order"):
            bwa._to_result(query=query, rec=rec)


def test_to_result_num_hits_on_unmapped(ref_fasta: Path) -> None:
    with BwaAlnInteractive(ref=ref_fasta, max_hits=1) as bwa:
        query = Query(bases="GATTACA", id="foo")
        # Exception: HN cannot be non-zero
        rec = SamBuilder().add_single(name=query.id, attrs={"HN": 42})
        with pytest.raises(ValueError, match="Read was unmapped but num_hits > 0"):
            bwa._to_result(query=query, rec=rec)
        # OK: HN tag is zero
        rec = SamBuilder().add_single(name=query.id, attrs={"HN": 0})
        bwa._to_result(query=query, rec=rec)
        # OK: no HN tag
        rec = SamBuilder().add_single(name=query.id)
        bwa._to_result(query=query, rec=rec)
