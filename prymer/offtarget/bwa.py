"""
# Methods and Classes to run and interact with BWA

The [`BwaAlnInteractive`][prymer.offtarget.bwa.BwaAlnInteractive] class is used to map
a list of queries to the reference genome.  A single query may be mapped using
[`map_one()`][prymer.offtarget.bwa.BwaAlnInteractive.map_one], while multiple queries may
be mapped using [`map_all()`][prymer.offtarget.bwa.BwaAlnInteractive.map_all].  The latter
is more efficient than mapping queries one-by-one.

The queries are provided to these methods using the [`Query`][prymer.offtarget.bwa.Query]
class, which contains the unique identifier for the query and the bases to query.  These methods
return a [`BwaResult`][prymer.offtarget.bwa.BwaResult], which represents zero or more hits
(or alignments) found by BWA for the given query.  Each hit is represented by a
[`BwaHit`][prymer.offtarget.bwa.BwaHit] object.  In some cases, BWA will write fewer (or no)
hits in the "XA" tag than the total number hits reported in the "HN".  This occurs when BWA finds more
hits than `max_hits` (see `bwt aln -X`).

 ## Example

```python
>>> from pathlib import Path
>>> ref_fasta = Path("./tests/offtarget/data/miniref.fa")
>>> query = Query(bases="TCTACTAAAAATACAAAAAATTAGCTGGGCATGATGGCATGCACCTGTAATCCCGCTACT", id="NA")
>>> bwa = BwaAlnInteractive(ref=ref_fasta, max_hits=1)
>>> result = bwa.map_one(query=query.bases, id=query.id)
>>> result.hit_count
1
>>> len(result.hits)
1
>>> result.hits[0]
BwaHit(refname='chr1', start=61, negative=False, cigar=Cigar(elements=(CigarElement(length=60, operator=<CigarOp.M: (0, 'M', True, True)>),)), edits=0)
>>> query = Query(bases="AAAAAA", id="NA")
>>> bwa.map_all(queries=[query])
[BwaResult(query=Query(id='NA', bases='AAAAAA'), hit_count=3968, hits=[])]
>>> bwa.close()

```
"""  # noqa: E501

import os
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar
from typing import Optional
from typing import cast

import pysam
from fgpyo import sam
from fgpyo import sequence
from fgpyo.sam import Cigar
from pysam import AlignedSegment

from prymer.api import coordmath
from prymer.util.executable_runner import ExecutableRunner


@dataclass(init=True, frozen=True)
class Query:
    """Represents a single query sequence for mapping.

    Attributes:
        id: the identifier for the query (e.g. read name)
        bases: the query bases
    """

    id: str
    bases: str

    _DEFAULT_BASE_QUALITY: ClassVar[str] = "H"
    """The default base quality for a query"""

    def __post_init__(self) -> None:
        if self.bases is None or len(self.bases) == 0:
            raise ValueError("No bases were provided to query")

    def to_fastq(self, reverse_complement: bool = False) -> str:
        """Returns the multi-line FASTQ string representation of this query"""
        bases = sequence.reverse_complement(self.bases) if reverse_complement else self.bases
        quals = self._DEFAULT_BASE_QUALITY * len(self.bases)
        return f"@{self.id}\n{bases}\n+\n{quals}\n"


@dataclass(init=True, frozen=True)
class BwaHit:
    """Represents a single hit or alignment of a sequence to a location in the genome.

    Attributes:
        refname: the reference name of the hit
        start: the start position of the hit (1-based inclusive)
        negative: whether the hit is on the negative strand
        cigar: the cigar string returned by BWA
        edits: the number of edits between the read and the reference
    """

    refname: str
    start: int
    negative: bool
    cigar: Cigar
    edits: int

    @property
    def mismatches(self) -> int:
        """The number of mismatches for the hit"""
        indel_sum = sum(elem.length for elem in self.cigar.elements if elem.operator.is_indel)
        if indel_sum > self.edits:
            raise ValueError(
                f"indel_sum ({indel_sum}) > self.edits ({self.edits}) with cigar: {self.cigar}"
            )
        return self.edits - indel_sum

    @property
    def end(self) -> int:
        """The end position of the hit (1-based inclusive)"""
        return coordmath.get_closed_end(self.start, self.cigar.length_on_target())

    @staticmethod
    def build(
        refname: str, start: int, negative: bool, cigar: str, edits: int, revcomp: bool = False
    ) -> "BwaHit":
        """Generates a hit object.

        Args:
            refname: the reference name of the hit
            start: the start position of the hit (1-based inclusive)
            negative: whether the hit is on the negative strand
            cigar: the cigar string returned by BWA
            edits: the number of edits between the read and the reference
            revcomp: whether the reverse-complement of the query sequence was fed to BWA, in which
                case various fields should be negated/reversed in the Hit

        Returns:
            A Hit that represents the mapping of the original query sequence that was supplied
        """
        negative = negative != revcomp
        return BwaHit(
            refname=refname,
            start=start,
            negative=negative,
            cigar=Cigar.from_cigarstring(cigar),
            edits=edits,
        )

    def __str__(self) -> str:
        """Returns the string representation in bwa's XA tag format."""
        # E.g. XA:Z:chr4,-97592047,24M,3;chr8,-32368749,24M,3;
        return ",".join(
            [
                self.refname,
                ("-" if self.negative else "+") + f"{self.start}",
                f"{self.cigar}",
                f"{self.edits}",
            ]
        )


@dataclass(init=True, frozen=True)
class BwaResult:
    """Represents zero or more hits or alignments found by BWA for the given query.  The number of
    hits found may be more than the number of hits listed in the `hits` attribute.

    Attributes:
        query: the query (as given, no RC'ing)
        hit_count: the total number of hits found by bwa (this may be more than len(hits))
        hits: the subset of hits that were reported by bwa
    """

    query: Query
    hit_count: int
    hits: list[BwaHit]


MAX_MISMATCHES: int = 3
"""The default maximum number of mismatches allowed in the full query sequence"""

MAX_MISMATCHES_IN_SEED: int = 3
"""The default maximum number of mismatches allowed in the seed region"""

MAX_GAP_OPENS: int = 0
"""The default maximum number of gap opens allowed in the full query sequence"""

MAX_GAP_EXTENSIONS: int = -1
"""The default maximum number of gap extensions allowed in the full query sequence"""

SEED_LENGTH: int = 20
"""The default length of the seed region"""

BWA_AUX_EXTENSIONS: list[str] = [".amb", ".ann", ".bwt", ".pac", ".sa"]
"""The file extensiosn for BWA index files"""


class BwaAlnInteractive(ExecutableRunner):
    """Wrapper around a novel mode of 'bwa aln' that allows for "interactive" use of bwa to keep
    the process running and be able to send it chunks of reads periodically and get alignments
    back without waiting for a full batch of reads to be sent.

    See: https://github.com/fulcrumgenomics/bwa/tree/interactive_aln

    Attributes:
        max_hits: the maximum number of hits to report - if more than this number of seed hits
                  are found, report only the count and not each hit.
        reverse_complement: reverse complement each query sequence before alignment.
        include_alt_hits: if True include hits to references with names ending in _alt, otherwise
                          do not include them.
    """

    def __init__(
        self,
        ref: Path,
        max_hits: int,
        executable: str | Path = "bwa",
        max_mismatches: int = 3,
        max_mismatches_in_seed: int = 3,
        max_gap_opens: int = 0,
        max_gap_extensions: int = -1,
        seed_length: int = 20,
        reverse_complement: bool = False,
        include_alt_hits: bool = False,
        threads: Optional[int] = None,
    ):
        """
        Args:
            ref: the path to the reference FASTA, which must be indexed with bwa.
            max_hits: the maximum number of hits to report - if more than this number of seed hits
                      are found, report only the count and not each hit.
            executable: string or Path representation of the `bwa` executable path
            max_mismatches: the maximum number of mismatches allowed in the full query sequence
            max_mismatches_in_seed: the maximum number of mismatches allowed in the seed region
            max_gap_opens: the maximum number of gap opens allowed in the full query sequence
            max_gap_extensions: the maximum number of gap extensions allowed in the full query
                                sequence
            seed_length: the length of the seed region
            reverse_complement: reverse complement each query sequence before alignment
            include_alt_hits: if true include hits to references with names ending in _alt,
                              otherwise do not include them.
            threads: the number of threads to use.  If `None`, use all available CPUs.
        """
        threads = os.cpu_count() if threads is None else threads
        executable_path = ExecutableRunner.validate_executable_path(executable=executable)
        self.reverse_complement: bool = reverse_complement
        self.include_alt_hits: bool = include_alt_hits
        self.max_hits: int = max_hits

        missing_aux_paths = []
        for aux_ext in BWA_AUX_EXTENSIONS:
            aux_path = Path(f"{ref}{aux_ext}")
            if not aux_path.exists():
                missing_aux_paths.append(aux_path)
        if len(missing_aux_paths) > 0:
            message: str
            if len(missing_aux_paths) > 1:
                message = "BWA index files do not exist:\n\t"
            else:
                message = "BWA index file does not exist:\n\t"
            message += "\t\n".join(f"{p}" for p in missing_aux_paths)
            raise FileNotFoundError(f"{message}\nPlease index with: `bwa index {ref}`")

        # -N = non-iterative mode: search for all n-difference hits (slooow)
        # -S = output SAM (run samse)
        # -Z = interactive mode (no input buffer and force processing with empty lines between recs)
        command: list[str] = [
            f"{executable_path}",
            "aln",
            "-t",
            f"{threads}",
            "-n",
            f"{max_mismatches}",
            "-o",
            f"{max_gap_opens}",
            "-e",
            f"{max_gap_extensions}",
            "-l",
            f"{seed_length}",
            "-k",
            f"{max_mismatches_in_seed}",
            "-X",
            f"{max_hits}",
            "-N",
            "-S",
            "-Z",
            "-D",
            f"{ref}",
            "/dev/stdin",
        ]

        super().__init__(command=command)

        # HACK ALERT
        # This is a hack. By trial and error, pysam.AlignmentFile() will block reading unless
        # there's at least a few records due to htslib wanting to read a few records for format
        # auto-detection. Lame.  So a hundred queries are sent to the aligner to align enable the
        # htslib auto-detection to complete, and for us to be able to read using pysam.
        num_warmup: int = 100
        for i in range(num_warmup):
            query = Query(id=f"ignoreme:{i}", bases="A" * 100)
            fastq_str = query.to_fastq(reverse_complement=self.reverse_complement)
            self._subprocess.stdin.write(fastq_str)
        self.__signal_bwa()  # forces the input to be sent to the underlying process.
        self._reader = sam.reader(path=self._subprocess.stdout, file_type=sam.SamFileType.SAM)
        for _ in range(num_warmup):
            next(self._reader)

    def __signal_bwa(self) -> None:
        """Signals BWA to process the queries"""
        for _ in range(3):
            self._subprocess.stdin.flush()
            self._subprocess.stdin.write("\n\n")
            self._subprocess.stdin.flush()

    def map_one(self, query: str, id: str = "unknown") -> BwaResult:
        """Maps a single query to the genome and returns the result.

        Args:
            query: the query to map with BWA

        Returns:
            a `BwaResult` based on mapping the query
        """
        return self.map_all([Query(bases=query, id=id)])[0]

    def map_all(self, queries: list[Query]) -> list[BwaResult]:
        """Maps multiple queries and returns the results.  This is more efficient than using
        `map_one` on each query one-by-one as it batches reads to BWA.

        Args:
            queries: the queries to map with BWA

        Returns:
            one `BwaResult`s for each query
        """
        if len(queries) == 0:
            return []

        # Send the reads to BWA
        for query in queries:
            fastq_str = query.to_fastq(reverse_complement=self.reverse_complement)
            self._subprocess.stdin.write(fastq_str)
        self.__signal_bwa()  # forces the input to be sent to the underlying process.

        # Read back the results
        results: list[BwaResult] = []
        for query in queries:
            # get the next alignment and convert to a result
            results.append(self._to_result(query=query, rec=next(self._reader)))

        return results

    def _to_result(self, query: Query, rec: pysam.AlignedSegment) -> BwaResult:
        """Converts the query and alignment to a result.

        Args:
            query: the original query
            rec: the alignment
        """
        if query.id != rec.query_name:
            raise ValueError(
                "Query and Results are out of order" f"Query=${query.id}, Result=${rec.query_name}"
            )

        num_hits: Optional[int] = int(rec.get_tag("HN")) if rec.has_tag("HN") else None
        if rec.is_unmapped:
            if num_hits is not None and num_hits > 0:
                raise ValueError(f"Read was unmapped but num_hits > 0: {rec}")
            return BwaResult(query=query, hit_count=0, hits=[])
        elif num_hits > self.max_hits:
            return BwaResult(query=query, hit_count=num_hits, hits=[])
        else:
            hits = self.to_hits(rec=rec)
            hit_count = num_hits if len(hits) == 0 else len(hits)
            return BwaResult(query=query, hit_count=hit_count, hits=hits)

    def to_hits(self, rec: AlignedSegment) -> list[BwaHit]:
        """Extracts the hits from the given alignment.  Beyond the current alignment
        additional alignments are parsed from the XA SAM tag.

        Args:
            rec: the given alignment
        """
        negative = rec.is_reverse != self.reverse_complement
        first_hit: BwaHit = BwaHit(
            refname=rec.reference_name,
            start=rec.reference_start + 1,  # NB: pysam is 0-based, Hit is 1-based
            negative=negative,
            cigar=Cigar.from_cigartuples(rec.cigartuples),
            edits=int(rec.get_tag("NM")),
        )

        # Add the hits
        # E.g. XA:Z:chr4,-97592047,24M,3;chr8,-32368749,24M,3;
        hits: list[BwaHit] = [first_hit]
        if rec.has_tag("XA"):
            for xa in cast(str, rec.get_tag("XA")).split(";"):
                if xa == "":
                    continue
                fields = xa.split(",")

                # If the reverse-complement of the query sequence was fed to BWA, various fields
                # should be negated/reversed in the Hit
                negative = fields[1][0] == "-"
                if self.reverse_complement:
                    negative = not negative

                hit: BwaHit = BwaHit(
                    refname=fields[0],
                    start=int(fields[1][1:]),  # the XA tag is 1-based
                    negative=negative,
                    cigar=Cigar.from_cigarstring(fields[2]),
                    edits=int(fields[3]),
                )
                hits.append(hit)

        if not self.include_alt_hits:
            hits = [hit for hit in hits if not hit.refname.endswith("_alt")]

        return hits

    def close(self) -> None:
        self._reader.close()
        super().close()
