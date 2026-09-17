#!/usr/bin/env python3

"""
Record reads whose BWA start positions hit annotated genomic features.
"""

from __future__ import annotations

import csv
import heapq
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Iterable


SCRIPT_DIR = Path(__file__).resolve().parent
ANNOTATION_DIR = SCRIPT_DIR / "gene_annotation"


@dataclass(frozen=True)
class ReadPosition:
    support_seq: str
    chrom: str
    pos1: int


Interval = tuple[int, int, str]


def open_text(path: Path):
    try:
        return path.open(encoding="utf-8", newline="")
    except OSError as exc:
        raise ValueError(f"cannot open {path}: {exc}") from exc


def load_reads(path: Path) -> list[ReadPosition]:
    reads: list[ReadPosition] = []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_id", "target_chrom", "insertion_pos", "strand"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "input header must contain sequence_id, target_chrom, "
                "insertion_pos and strand columns"
            )
        for line_number, row in enumerate(reader, 2):
            try:
                pos1 = int(row["insertion_pos"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: insertion_pos must be an integer"
                ) from exc
            if pos1 < 1:
                continue

            strand = row["strand"]
            if strand not in {"positive", "negative"}:
                raise ValueError(
                    f"{path}:{line_number}: strand must be positive or negative"
                )
            query_pos1 = pos1 if strand == "positive" else pos1 + 1
            reads.append(
                ReadPosition(row["sequence_id"], row["target_chrom"], query_pos1)
            )
    return reads


def bed_records(path: Path) -> Iterable[tuple[str, int, int, list[str]]]:
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least 4 columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: BED start/end must be integers"
                ) from exc
            if start < 0 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            yield fields[0], start, end, fields[3:]


def load_gene_intervals(filename: str) -> dict[str, list[Interval]]:
    intervals: DefaultDict[str, list[Interval]] = defaultdict(list)
    path = ANNOTATION_DIR / filename
    for chrom, start, end, extra in bed_records(path):
        intervals[chrom].append((start, end, extra[0]))
    for chrom_intervals in intervals.values():
        chrom_intervals.sort()
    return dict(intervals)


def load_tss_gene_map() -> dict[tuple[str, str, str, int], set[str]]:
    path = ANNOTATION_DIR / "TSS_transcript_level.tsv"
    mapping: DefaultDict[tuple[str, str, str, int], set[str]] = defaultdict(set)
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "gene_id", "chrom", "strand", "TSS_pos"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path}: header must contain transcript_id, gene_id, chrom, "
                "strand and TSS_pos"
            )
        for line_number, row in enumerate(reader, 2):
            try:
                tss0 = int(row["TSS_pos"]) - 1
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{line_number}: TSS_pos must be an integer") from exc
            key = (row["transcript_id"], row["chrom"], row["strand"], tss0)
            mapping[key].add(row["gene_id"])
    return dict(mapping)


def load_upstream_intervals() -> dict[str, list[Interval]]:
    tss_to_gene = load_tss_gene_map()
    path = ANNOTATION_DIR / "TSS_transcript_level.0based.bed"
    intervals: DefaultDict[str, list[Interval]] = defaultdict(list)
    for chrom, start, end, extra in bed_records(path):
        if len(extra) < 3 or extra[2] not in {"+", "-"}:
            raise ValueError(f"{path}: expected strand (+/-) in BED column 6")
        transcript_id, strand = extra[0], extra[2]
        key = (transcript_id, chrom, strand, start)
        try:
            gene_ids = tss_to_gene[key]
        except KeyError as exc:
            raise ValueError(
                f"{path}: no gene_id mapping for TSS record {key!r}"
            ) from exc
        if strand == "+":
            upstream_start, upstream_end = max(0, start - 2000), start
        else:
            upstream_start, upstream_end = end, end + 2000
        if upstream_start < upstream_end:
            for gene_id in gene_ids:
                intervals[chrom].append((upstream_start, upstream_end, gene_id))
    for chrom_intervals in intervals.values():
        chrom_intervals.sort()
    return dict(intervals)


def find_hits(
    reads: list[ReadPosition], intervals: dict[str, list[Interval]]
) -> dict[int, set[str]]:
    """Find all genes overlapping each read using a chromosome-wise sweep."""
    queries: DefaultDict[str, list[tuple[int, int]]] = defaultdict(list)
    for read_number, read in enumerate(reads):
        queries[read.chrom].append((read.pos1 - 1, read_number))

    hits: dict[int, set[str]] = {}
    for chrom, chrom_queries in queries.items():
        chrom_intervals = intervals.get(chrom)
        if not chrom_intervals:
            continue
        active_ends: list[tuple[int, int, str]] = []
        active_genes: Counter[str] = Counter()
        interval_number = 0
        for pos0, read_number in sorted(chrom_queries):
            while (
                interval_number < len(chrom_intervals)
                and chrom_intervals[interval_number][0] <= pos0
            ):
                _, end, gene_id = chrom_intervals[interval_number]
                heapq.heappush(active_ends, (end, interval_number, gene_id))
                active_genes[gene_id] += 1
                interval_number += 1
            while active_ends and active_ends[0][0] <= pos0:
                _, _, gene_id = heapq.heappop(active_ends)
                active_genes[gene_id] -= 1
                if active_genes[gene_id] == 0:
                    del active_genes[gene_id]
            if active_genes:
                hits[read_number] = set(active_genes)
    return hits


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: python {Path(sys.argv[0]).name} [input]", file=sys.stderr)
        return 2
    try:
        reads = load_reads(Path(sys.argv[1]))
        exon_hits = find_hits(
            reads, load_gene_intervals("gene_exon_gene_level.0based.bed")
        )
        intron_hits = find_hits(
            reads, load_gene_intervals("gene_intron_gene_level.0based.bed")
        )
        upstream_hits = find_hits(reads, load_upstream_intervals())
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("gene_id", "chrom", "type", "insertion_site", "support_seq"))
    for read_number, read in enumerate(reads):
        if read_number in exon_hits:
            feature, gene_ids = "Exon", exon_hits[read_number]
        elif read_number in intron_hits:
            feature, gene_ids = "Intron", intron_hits[read_number]
        elif read_number in upstream_hits:
            feature, gene_ids = "Upstream TSS (2kb)", upstream_hits[read_number]
        else:
            continue
        writer.writerow(
            (",".join(sorted(gene_ids)), read.chrom, feature, read.pos1, read.support_seq)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
