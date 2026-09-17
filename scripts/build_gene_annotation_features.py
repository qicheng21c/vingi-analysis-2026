#!/usr/bin/env python3
"""
Build reusable hg38 gene annotation feature files from UCSC ncbiRefSeq-like tables.

Input coordinate assumption for UCSC ncbiRefSeq:
  - txStart / exonStarts: 0-based starts
  - txEnd / exonEnds: 0-based half-open ends

Output:
  - TSV feature coordinates: 1-based inclusive
  - BED coordinates: 0-based half-open
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


@dataclass
class Transcript:
    gene_id: str
    gene_name: str
    transcript_id: str
    chrom: str
    strand: str
    tx_start: int  # 1-based inclusive
    tx_end: int    # 1-based inclusive
    exons: List[Tuple[int, int]]  # 1-based inclusive
    warning: str = ""


@dataclass
class GeneFeature:
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    gene_start: int
    gene_end: int
    transcript_ids: List[str] = field(default_factory=list)
    exons: List[Tuple[int, int]] = field(default_factory=list)
    introns: List[Tuple[int, int]] = field(default_factory=list)
    has_exon_annotation: bool = True
    warning: str = ""

    @property
    def gene_length(self) -> int:
        return self.gene_end - self.gene_start + 1


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def ensure_output_dir(path: str) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def check_output_paths(paths: Sequence[Path], force: bool) -> None:
    existing = [str(p) for p in paths if p.exists()]
    if existing and not force:
        msg = "Output file(s) already exist. Re-run with --force to overwrite:\n" + "\n".join(existing)
        raise FileExistsError(msg)


def normalize_missing(value: str, fallback: str) -> str:
    value = (value or "").strip()
    if value in {"", ".", "NA", "N/A", "None", "none"}:
        return fallback
    return value


def safe_gene_id(gene_name: str, chrom: str, strand: str) -> str:
    base = normalize_missing(gene_name, "unknown_gene")
    return f"{base}|{chrom}|{strand}"


def parse_int(value: str, field_name: str, line_no: int) -> Optional[int]:
    try:
        return int(value)
    except Exception:
        sys.stderr.write(f"WARNING line {line_no}: non-integer {field_name}={value!r}; record skipped\n")
        return None


def split_comma_ints(value: str, line_no: int, field_name: str) -> Optional[List[int]]:
    value = (value or "").strip().rstrip(",")
    if value == "":
        return []
    out: List[int] = []
    for item in value.split(","):
        item = item.strip()
        if item == "":
            continue
        try:
            out.append(int(item))
        except Exception:
            sys.stderr.write(f"WARNING line {line_no}: non-integer in {field_name}: {item!r}; exon structure skipped\n")
            return None
    return out


def is_header(fields: Sequence[str]) -> bool:
    lowered = {x.strip().lstrip("#") for x in fields}
    required = {"chrom", "strand"}
    return bool(required & lowered) and ("txStart" in lowered or "txEnd" in lowered or "exonStarts" in lowered)


def make_index_from_header(fields: Sequence[str]) -> Dict[str, int]:
    idx = {}
    for i, f in enumerate(fields):
        name = f.strip().lstrip("#")
        idx[name] = i
    return idx


def make_index_no_header(n_fields: int) -> Dict[str, int]:
    """Return column indices for UCSC ncbiRefSeq/refGene-like rows with or without a leading bin column."""
    # With bin: bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
    # exonCount, exonStarts, exonEnds, score, name2, ...
    # Without bin: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
    # exonCount, exonStarts, exonEnds, score, name2, ...
    if n_fields >= 13 and re.fullmatch(r"\d+", "0"):
        # Decide from the third column looking like chromosome and fourth looking like strand.
        # This function does not see the row values, so caller may override by passing row-aware index.
        pass
    raise RuntimeError("make_index_no_header should be called through make_index_from_row")


def make_index_from_row(fields: Sequence[str]) -> Dict[str, int]:
    if len(fields) >= 11 and len(fields[3]) == 1 and fields[3] in {"+", "-"}:
        # leading bin column present
        return {
            "name": 1,
            "chrom": 2,
            "strand": 3,
            "txStart": 4,
            "txEnd": 5,
            "exonCount": 8,
            "exonStarts": 9,
            "exonEnds": 10,
            "name2": 12 if len(fields) > 12 else 1,
        }
    if len(fields) >= 10 and len(fields[2]) == 1 and fields[2] in {"+", "-"}:
        # no bin column
        return {
            "name": 0,
            "chrom": 1,
            "strand": 2,
            "txStart": 3,
            "txEnd": 4,
            "exonCount": 7,
            "exonStarts": 8,
            "exonEnds": 9,
            "name2": 11 if len(fields) > 11 else 0,
        }
    raise ValueError(
        "Cannot infer UCSC ncbiRefSeq/refGene columns. Expected either a row with leading bin "
        "or no leading bin."
    )


def get_field(fields: Sequence[str], index: Dict[str, int], name: str, default: str = "") -> str:
    i = index.get(name)
    if i is None or i >= len(fields):
        return default
    return fields[i]


def parse_gene_annotation(path: str) -> Tuple[List[Transcript], Dict[str, int]]:
    transcripts: List[Transcript] = []
    report = defaultdict(int)
    header_index: Optional[Dict[str, int]] = None
    inferred_index: Optional[Dict[str, int]] = None

    with open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                report["empty_lines"] += 1
                continue
            if line.startswith("#") and not is_header(line.split("\t")):
                report["comment_lines"] += 1
                continue

            fields = line.split("\t")
            if header_index is None and inferred_index is None and is_header(fields):
                header_index = make_index_from_header(fields)
                report["header_lines"] += 1
                continue

            try:
                index = header_index if header_index is not None else inferred_index
                if index is None:
                    index = make_index_from_row(fields)
                    inferred_index = index
            except Exception as exc:
                report["skipped_unparseable_column_layout"] += 1
                sys.stderr.write(f"WARNING line {line_no}: {exc}; record skipped\n")
                continue

            transcript_id = normalize_missing(get_field(fields, index, "name"), f"transcript_line_{line_no}")
            chrom = normalize_missing(get_field(fields, index, "chrom"), "unknown_chrom")
            strand = normalize_missing(get_field(fields, index, "strand"), ".")
            gene_name = normalize_missing(get_field(fields, index, "name2", transcript_id), transcript_id)

            if strand not in {"+", "-"}:
                report["skipped_missing_or_invalid_strand"] += 1
                sys.stderr.write(f"WARNING line {line_no}: invalid strand={strand!r}; record skipped\n")
                continue

            tx_start0 = parse_int(get_field(fields, index, "txStart"), "txStart", line_no)
            tx_end0 = parse_int(get_field(fields, index, "txEnd"), "txEnd", line_no)
            if tx_start0 is None or tx_end0 is None or tx_end0 <= tx_start0:
                report["skipped_invalid_transcript_interval"] += 1
                continue

            tx_start1 = tx_start0 + 1
            tx_end1 = tx_end0

            exon_starts0 = split_comma_ints(get_field(fields, index, "exonStarts"), line_no, "exonStarts")
            exon_ends0 = split_comma_ints(get_field(fields, index, "exonEnds"), line_no, "exonEnds")
            exon_count_raw = parse_int(get_field(fields, index, "exonCount", "0"), "exonCount", line_no)

            warnings: List[str] = []
            exons: List[Tuple[int, int]] = []
            if exon_starts0 is None or exon_ends0 is None:
                warnings.append("invalid_exon_structure")
                report["transcripts_with_invalid_exon_structure"] += 1
            elif len(exon_starts0) != len(exon_ends0):
                warnings.append("exon_start_end_count_mismatch")
                report["transcripts_with_exon_start_end_count_mismatch"] += 1
            else:
                if exon_count_raw is not None and exon_count_raw != len(exon_starts0):
                    warnings.append("exonCount_mismatch")
                    report["transcripts_with_exonCount_mismatch"] += 1
                for s0, e0 in zip(exon_starts0, exon_ends0):
                    if e0 <= s0:
                        warnings.append("invalid_exon_interval")
                        report["invalid_exon_intervals"] += 1
                        continue
                    s1 = max(tx_start1, s0 + 1)
                    e1 = min(tx_end1, e0)
                    if e1 >= s1:
                        exons.append((s1, e1))
                if not exons:
                    warnings.append("no_valid_exon_annotation")
                    report["transcripts_without_valid_exon_annotation"] += 1

            transcripts.append(
                Transcript(
                    gene_id=safe_gene_id(gene_name, chrom, strand),
                    gene_name=gene_name,
                    transcript_id=transcript_id,
                    chrom=chrom,
                    strand=strand,
                    tx_start=tx_start1,
                    tx_end=tx_end1,
                    exons=exons,
                    warning=";".join(sorted(set(warnings))),
                )
            )
            report["parsed_transcripts"] += 1

    report["total_transcripts_returned"] = len(transcripts)
    return transcripts, dict(report)


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    sorted_intervals = sorted((s, e) for s, e in intervals if e >= s)
    if not sorted_intervals:
        return []
    merged: List[Tuple[int, int]] = []
    cur_s, cur_e = sorted_intervals[0]
    for s, e in sorted_intervals[1:]:
        if s <= cur_e + 1:  # merge overlapping or adjacent intervals
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def subtract_intervals(body: Tuple[int, int], exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    body_start, body_end = body
    introns: List[Tuple[int, int]] = []
    cursor = body_start
    for s, e in sorted(exons):
        if s > cursor:
            introns.append((cursor, s - 1))
        cursor = max(cursor, e + 1)
    if cursor <= body_end:
        introns.append((cursor, body_end))
    return [(s, e) for s, e in introns if e >= s]


def build_gene_features(transcripts: Sequence[Transcript]) -> List[GeneFeature]:
    grouped: Dict[str, List[Transcript]] = defaultdict(list)
    for tr in transcripts:
        grouped[tr.gene_id].append(tr)

    genes: List[GeneFeature] = []
    for gene_id, trs in grouped.items():
        trs_sorted = sorted(trs, key=lambda x: (x.chrom, x.tx_start, x.tx_end, x.transcript_id))
        first = trs_sorted[0]
        gene_start = min(t.tx_start for t in trs_sorted)
        gene_end = max(t.tx_end for t in trs_sorted)
        warnings = sorted({t.warning for t in trs_sorted if t.warning})
        all_exons = []
        for tr in trs_sorted:
            all_exons.extend(tr.exons)
        merged_exons = merge_intervals(all_exons)
        has_exon = bool(merged_exons)
        if not has_exon:
            warnings.append("gene_has_no_valid_exon_annotation")
        introns = subtract_intervals((gene_start, gene_end), merged_exons) if has_exon else []
        genes.append(
            GeneFeature(
                gene_id=gene_id,
                gene_name=first.gene_name,
                chrom=first.chrom,
                strand=first.strand,
                gene_start=gene_start,
                gene_end=gene_end,
                transcript_ids=[t.transcript_id for t in trs_sorted],
                exons=merged_exons,
                introns=introns,
                has_exon_annotation=has_exon,
                warning=";".join(sorted(set(warnings))),
            )
        )
    return sorted(genes, key=lambda g: (g.chrom, g.gene_start, g.gene_end, g.gene_name, g.strand))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> int:
    n = 0
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            n += 1
    return n


def write_raw_transcripts(path: Path, transcripts: Sequence[Transcript], source: str, genome: str) -> int:
    fields = [
        "gene_id", "gene_name", "transcript_id", "chrom", "strand",
        "transcript_start", "transcript_end", "exon_count",
        "exon_starts_1based", "exon_ends_1based", "source_annotation",
        "genome_build", "warning",
    ]
    rows = []
    for tr in transcripts:
        rows.append({
            "gene_id": tr.gene_id,
            "gene_name": tr.gene_name,
            "transcript_id": tr.transcript_id,
            "chrom": tr.chrom,
            "strand": tr.strand,
            "transcript_start": tr.tx_start,
            "transcript_end": tr.tx_end,
            "exon_count": len(tr.exons),
            "exon_starts_1based": ",".join(str(s) for s, _ in tr.exons),
            "exon_ends_1based": ",".join(str(e) for _, e in tr.exons),
            "source_annotation": source,
            "genome_build": genome,
            "warning": tr.warning,
        })
    return write_tsv(path, fields, rows)


def write_gene_body(path: Path, genes: Sequence[GeneFeature], source: str, genome: str) -> int:
    fields = [
        "gene_id", "gene_name", "chrom", "strand", "gene_start", "gene_end",
        "gene_length", "transcript_count", "source_annotation", "genome_build",
    ]
    rows = []
    for g in genes:
        rows.append({
            "gene_id": g.gene_id,
            "gene_name": g.gene_name,
            "chrom": g.chrom,
            "strand": g.strand,
            "gene_start": g.gene_start,
            "gene_end": g.gene_end,
            "gene_length": g.gene_length,
            "transcript_count": len(g.transcript_ids),
            "source_annotation": source,
            "genome_build": genome,
        })
    return write_tsv(path, fields, rows)


def write_feature_intervals(path: Path, genes: Sequence[GeneFeature]) -> int:
    fields = [
        "gene_id", "gene_name", "chrom", "strand", "gene_start", "gene_end",
        "feature_type", "feature_index", "feature_start", "feature_end", "feature_length",
    ]
    rows = []
    for g in genes:
        for feature_type, intervals in (("Exon", g.exons), ("Intron", g.introns)):
            for i, (s, e) in enumerate(intervals, start=1):
                rows.append({
                    "gene_id": g.gene_id,
                    "gene_name": g.gene_name,
                    "chrom": g.chrom,
                    "strand": g.strand,
                    "gene_start": g.gene_start,
                    "gene_end": g.gene_end,
                    "feature_type": feature_type,
                    "feature_index": i,
                    "feature_start": s,
                    "feature_end": e,
                    "feature_length": e - s + 1,
                })
    return write_tsv(path, fields, rows)


def write_summary(path: Path, genes: Sequence[GeneFeature]) -> int:
    fields = [
        "gene_id", "gene_name", "chrom", "strand", "gene_start", "gene_end",
        "gene_length", "transcript_count", "exon_interval_count", "intron_interval_count",
        "total_exon_length", "total_intron_length", "exon_coverage_fraction",
        "intron_coverage_fraction", "has_exon_annotation", "warning",
    ]
    rows = []
    for g in genes:
        total_exon = sum(e - s + 1 for s, e in g.exons)
        total_intron = sum(e - s + 1 for s, e in g.introns)
        rows.append({
            "gene_id": g.gene_id,
            "gene_name": g.gene_name,
            "chrom": g.chrom,
            "strand": g.strand,
            "gene_start": g.gene_start,
            "gene_end": g.gene_end,
            "gene_length": g.gene_length,
            "transcript_count": len(g.transcript_ids),
            "exon_interval_count": len(g.exons),
            "intron_interval_count": len(g.introns),
            "total_exon_length": total_exon,
            "total_intron_length": total_intron,
            "exon_coverage_fraction": f"{total_exon / g.gene_length:.6f}" if g.gene_length else "0.000000",
            "intron_coverage_fraction": f"{total_intron / g.gene_length:.6f}" if g.gene_length else "0.000000",
            "has_exon_annotation": "true" if g.has_exon_annotation else "false",
            "warning": g.warning,
        })
    return write_tsv(path, fields, rows)


def write_tss(path: Path, transcripts: Sequence[Transcript], genome: str) -> int:
    fields = [
        "gene_id", "gene_name", "transcript_id", "chrom", "strand",
        "transcript_start", "transcript_end", "TSS_pos", "genome_build",
    ]
    rows = []
    for tr in transcripts:
        tss = tr.tx_start if tr.strand == "+" else tr.tx_end
        rows.append({
            "gene_id": tr.gene_id,
            "gene_name": tr.gene_name,
            "transcript_id": tr.transcript_id,
            "chrom": tr.chrom,
            "strand": tr.strand,
            "transcript_start": tr.tx_start,
            "transcript_end": tr.tx_end,
            "TSS_pos": tss,
            "genome_build": genome,
        })
    rows.sort(key=lambda r: (r["chrom"], int(r["TSS_pos"]), r["gene_name"], r["transcript_id"]))
    return write_tsv(path, fields, rows)


def write_bed(path: Path, rows: Iterable[Tuple[str, int, int, str, str, str]]) -> int:
    n = 0
    with open(path, "w", encoding="utf-8", newline="") as fh:
        for chrom, start1, end1, name, score, strand in rows:
            start0 = max(0, start1 - 1)
            end0 = end1
            fh.write(f"{chrom}\t{start0}\t{end0}\t{name}\t{score}\t{strand}\n")
            n += 1
    return n


def write_report(path: Path, report_rows: Sequence[Tuple[str, object]]) -> int:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("item\tvalue\n")
        for key, value in report_rows:
            fh.write(f"{key}\t{value}\n")
    return len(report_rows)


def build_outputs(args: argparse.Namespace) -> None:
    out = ensure_output_dir(args.output_dir)
    paths = {
        "raw": out / "gene_annotation.raw.tsv",
        "body": out / "gene_body_gene_level.tsv",
        "body_bed": out / "gene_body_gene_level.0based.bed",
        "feature": out / "gene_exon_intron_gene_level.tsv",
        "exon_bed": out / "gene_exon_gene_level.0based.bed",
        "intron_bed": out / "gene_intron_gene_level.0based.bed",
        "summary": out / "gene_exon_intron_summary.tsv",
        "tss": out / "TSS_transcript_level.tsv",
        "tss_bed": out / "TSS_transcript_level.0based.bed",
        "report": out / "annotation_build_report.tsv",
    }
    check_output_paths(list(paths.values()), args.force)

    transcripts, parse_report = parse_gene_annotation(args.gene_annotation)
    if not transcripts:
        raise RuntimeError("No valid transcript records were parsed from the gene annotation file.")

    genes = build_gene_features(transcripts)

    counts = {}
    counts["gene_annotation.raw.tsv_rows"] = write_raw_transcripts(paths["raw"], transcripts, args.source_annotation, args.genome_build)
    counts["gene_body_gene_level.tsv_rows"] = write_gene_body(paths["body"], genes, args.source_annotation, args.genome_build)
    counts["gene_body_gene_level.0based.bed_rows"] = write_bed(
        paths["body_bed"],
        ((g.chrom, g.gene_start, g.gene_end, g.gene_id, "0", g.strand) for g in genes),
    )
    counts["gene_exon_intron_gene_level.tsv_rows"] = write_feature_intervals(paths["feature"], genes)
    counts["gene_exon_gene_level.0based.bed_rows"] = write_bed(
        paths["exon_bed"],
        ((g.chrom, s, e, g.gene_id, "0", g.strand) for g in genes for s, e in g.exons),
    )
    counts["gene_intron_gene_level.0based.bed_rows"] = write_bed(
        paths["intron_bed"],
        ((g.chrom, s, e, g.gene_id, "0", g.strand) for g in genes for s, e in g.introns),
    )
    counts["gene_exon_intron_summary.tsv_rows"] = write_summary(paths["summary"], genes)
    counts["TSS_transcript_level.tsv_rows"] = write_tss(paths["tss"], transcripts, args.genome_build)
    counts["TSS_transcript_level.0based.bed_rows"] = write_bed(
        paths["tss_bed"],
        ((tr.chrom, tr.tx_start if tr.strand == "+" else tr.tx_end,
          tr.tx_start if tr.strand == "+" else tr.tx_end,
          tr.transcript_id, "0", tr.strand) for tr in transcripts),
    )

    report_rows: List[Tuple[str, object]] = [
        ("script", Path(sys.argv[0]).name),
        ("input_gene_annotation", args.gene_annotation),
        ("output_dir", args.output_dir),
        ("source_annotation", args.source_annotation),
        ("genome_build", args.genome_build),
        ("coordinate_system_tsv", "1-based inclusive"),
        ("coordinate_system_bed", "0-based half-open"),
        ("gene_count", len(genes)),
        ("transcript_count", len(transcripts)),
        ("genes_without_valid_exon_annotation", sum(1 for g in genes if not g.has_exon_annotation)),
    ]
    report_rows.extend((f"parse_{k}", v) for k, v in sorted(parse_report.items()))
    report_rows.extend((k, v) for k, v in sorted(counts.items()))
    write_report(paths["report"], report_rows)

    sys.stderr.write("Finished building annotation feature files.\n")
    for name, path in paths.items():
        sys.stderr.write(f"  {name}: {path}\n")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build gene-level body/exon/intron and transcript-level TSS files from UCSC ncbiRefSeq-like annotation."
    )
    parser.add_argument("--gene-annotation", required=True, help="Input gene annotation file, e.g. UCSC ncbiRefSeq.txt.gz")
    parser.add_argument("--output-dir", required=True, help="Output directory for reusable annotation feature files")
    parser.add_argument("--source-annotation", default="UCSC ncbiRefSeq", help="Source annotation label written to output tables")
    parser.add_argument("--genome-build", default="hg38", help="Genome build label written to output tables")
    parser.add_argument("--force", action="store_true", help="Overwrite output files created by this script")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        build_outputs(args)
    except Exception as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

