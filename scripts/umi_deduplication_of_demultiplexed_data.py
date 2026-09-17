#!/usr/bin/env python3

import argparse
import gzip
import math
import os
import re
import sys
from contextlib import ExitStack
from dataclasses import dataclass
from typing import TextIO


UMI_PATTERN = re.compile(r"(?:^|\s)umi=([^\s]+)")


@dataclass(frozen=True)
class FastqRecord:
    """Store FASTQ content read from up to four lines."""

    lines: tuple[str, ...]

    @property
    def header(self) -> str:
        return self.lines[0] if self.lines else ""

    @property
    def sequence(self) -> str:
        return self.lines[1] if len(self.lines) > 1 else ""

    @property
    def separator(self) -> str:
        return self.lines[2] if len(self.lines) > 2 else ""

    @property
    def quality(self) -> str:
        return self.lines[3] if len(self.lines) > 3 else ""


@dataclass
class BestPair:
    """Store the best read pair for a given combination of UMI and 6 bp after the adapter."""

    score: float
    input_order: int
    r1: FastqRecord
    r2: FastqRecord


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Deduplicate paired-end FASTQ files by UMI and 6 bp after the adapter, automatically skipping invalid records."
    )
    parser.add_argument("input_r1", help="R1 FASTQ input file")
    parser.add_argument("input_r2", help="R2 FASTQ input file")
    parser.add_argument("output_r1", help="Deduplicated R1 FASTQ output file")
    parser.add_argument("output_r2", help="Deduplicated R2 FASTQ output file")
    return parser.parse_args()


def open_text(path: str, mode: str) -> TextIO:
    """Open a plain-text or gzip-compressed file based on its extension."""

    if path.lower().endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def read_record(handle: TextIO) -> FastqRecord | None:
    """Read up to four lines from the current position and return None at end of file."""

    first_line = handle.readline()
    if first_line == "":
        return None

    lines = [first_line]
    for _ in range(3):
        line = handle.readline()
        if line == "":
            break
        lines.append(line)
    return FastqRecord(tuple(lines))


def extract_umi(header: str) -> str | None:
    """Extract a non-empty UMI from the sequence header."""

    match = UMI_PATTERN.search(header.rstrip("\r\n"))
    return match.group(1) if match is not None else None


def read_id(header: str) -> str:
    """Extract the instrument read ID before whitespace, ignoring the leading '@'."""

    return header[1:].split(maxsplit=1)[0]


def validate_record(record: FastqRecord, side: str) -> list[str]:
    """Perform six format checks on a single-end FASTQ record."""

    errors: list[str] = []
    if len(record.lines) != 4:
        errors.append(f"{side} record is incomplete, with {len(record.lines)} lines")
        return errors
    if not record.header.startswith("@"):
        errors.append(f"{side} first line does not start with '@'")
    if not record.separator.startswith("+"):
        errors.append(f"{side} third line does not start with '+'")

    sequence = record.sequence.rstrip("\r\n")
    quality = record.quality.rstrip("\r\n")
    if len(sequence) != len(quality):
        errors.append(
            f"{side} sequence length {len(sequence)} does not match quality length {len(quality)}"
        )
    if not quality:
        errors.append(f"{side} quality line is empty")
    elif any(ord(character) < 33 or ord(character) > 126 for character in quality):
        errors.append(f"{side} contains quality characters outside the Phred+33 range")
    if extract_umi(record.header) is None:
        errors.append(f"{side} sequence header does not contain a non-empty umi=... field")
    return errors


def validate_pair(
    r1: FastqRecord | None, r2: FastqRecord | None
) -> tuple[list[str], str | None]:
    """Validate a read pair and return the UMI if valid."""

    errors: list[str] = []
    if r1 is None or r2 is None:
        errors.append("R1 and R2 contain different numbers of records")
    if r1 is not None:
        errors.extend(validate_record(r1, "R1"))
    if r2 is not None:
        errors.extend(validate_record(r2, "R2"))
    if errors or r1 is None or r2 is None:
        return errors, None

    r1_sequence = r1.sequence.rstrip("\r\n")
    if len(r1_sequence) < 6:
        errors.append(f"R1 sequence length {len(r1_sequence)} is less than 6 bp")
        return errors, None

    if read_id(r1.header) != read_id(r2.header):
        errors.append(
            f"R1/R2 read IDs do not match: {read_id(r1.header)} != {read_id(r2.header)}"
        )
    umi_r1 = extract_umi(r1.header)
    umi_r2 = extract_umi(r2.header)
    if umi_r1 != umi_r2:
        errors.append(f"R1/R2 UMIs do not match: {umi_r1} != {umi_r2}")
    return errors, umi_r1 if not errors else None


def mean_quality(record: FastqRecord) -> float:
    """Calculate the mean Phred+33 quality from the fourth line of a validated record."""

    quality = record.quality.rstrip("\r\n")
    return math.fsum(ord(character) - 33 for character in quality) / len(quality)


def validate_paths(args: argparse.Namespace) -> None:
    """Prevent output files from overwriting input files or each other."""

    input_paths = {os.path.abspath(args.input_r1), os.path.abspath(args.input_r2)}
    output_paths = {os.path.abspath(args.output_r1), os.path.abspath(args.output_r2)}
    if len(output_paths) != 2:
        raise ValueError("output_r1 and output_r2 must be different files")
    if input_paths & output_paths:
        raise ValueError("Output file paths must not match any input file path")


def select_best_pairs(
    args: argparse.Namespace,
) -> tuple[dict[tuple[str, str], BestPair], int, int]:
    """Scan input files and select the best read pair for each combination of UMI and 6 bp after the adapter."""

    best_by_molecule: dict[tuple[str, str], BestPair] = {}
    input_count = 0
    invalid_count = 0

    with ExitStack() as stack:
        r1_handle = stack.enter_context(open_text(args.input_r1, "rt"))
        r2_handle = stack.enter_context(open_text(args.input_r2, "rt"))
        while True:
            r1 = read_record(r1_handle)
            r2 = read_record(r2_handle)
            if r1 is None and r2 is None:
                break

            input_count += 1
            errors, umi = validate_pair(r1, r2)
            if errors:
                invalid_count += 1
                print(
                    f"Warning: skipping record {input_count}: {'; '.join(errors)}",
                    file=sys.stderr,
                )
                continue

            if r1 is None or r2 is None or umi is None:
                raise RuntimeError("Internal error: validated read pair is missing required data")
            adapter_after_6bp = r1.sequence.rstrip("\r\n")[:6]
            molecule_key = (umi, adapter_after_6bp)
            pair_score = (mean_quality(r1) + mean_quality(r2)) / 2.0
            current = best_by_molecule.get(molecule_key)
            if current is None or pair_score > current.score:
                best_by_molecule[molecule_key] = BestPair(
                    pair_score, input_count, r1, r2
                )

    return best_by_molecule, input_count, invalid_count


def format_output_header(record: FastqRecord, umi: str) -> str:
    """Generate the output header in the format @original_sequence_name:UMI:UMI_sequence."""

    return f"@{read_id(record.header)}:UMI:{umi}\n"


def write_record(handle: TextIO, record: FastqRecord, umi: str) -> None:
    """Write a four-line FASTQ record using the modified sequence header."""

    handle.write(format_output_header(record, umi))
    for line in record.lines[1:]:
        handle.write(line)


def write_best_pairs(
    args: argparse.Namespace,
    best_by_molecule: dict[tuple[str, str], BestPair],
) -> None:
    """Write read pairs according to the input positions of the final retained records."""

    selected = sorted(best_by_molecule.values(), key=lambda pair: pair.input_order)
    with ExitStack() as stack:
        r1_output = stack.enter_context(open_text(args.output_r1, "wt"))
        r2_output = stack.enter_context(open_text(args.output_r2, "wt"))
        for pair in selected:
            umi = extract_umi(pair.r1.header)
            if umi is None:
                raise RuntimeError("Internal error: output read pair is missing a UMI")
            write_record(r1_output, pair.r1, umi)
            write_record(r2_output, pair.r2, umi)


def main() -> int:
    args = parse_args()
    try:
        validate_paths(args)
        best_by_molecule, input_count, invalid_count = select_best_pairs(args)
        write_best_pairs(args, best_by_molecule)
    except (OSError, UnicodeError, EOFError, ValueError, RuntimeError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1

    valid_count = input_count - invalid_count
    duplicate_count = valid_count - len(best_by_molecule)
    print(
        f"Processing complete: read {input_count} record pairs, {valid_count} valid pairs, "
        f"skipped {invalid_count} invalid pairs, removed {duplicate_count} molecular duplicate pairs, "
        f"output {len(best_by_molecule)} read pairs.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())