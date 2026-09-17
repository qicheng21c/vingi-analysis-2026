#!/usr/bin/env python3
import csv
import re
import sys


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
MATCH_OPS = {"M", "=", "X"}


def matched_length_from_cigar(cigar, line_number):
    """Fully parse the CIGAR string and return the total length of M, =, and X operations."""
    if not isinstance(cigar, str) or not cigar:
        raise ValueError(f"cigar is empty at line {line_number}")

    ops = [(int(length), op) for length, op in CIGAR_PATTERN.findall(cigar)]
    if not ops or any(length <= 0 for length, _ in ops):
        raise ValueError(f"Invalid cigar at line {line_number}: {cigar!r}")

    reconstructed = "".join(f"{length}{op}" for length, op in ops)
    if reconstructed != cigar:
        raise ValueError(
            f"Incomplete cigar format or unknown operation at line {line_number}: {cigar!r}"
        )

    matched_length = sum(
        length for length, op in ops if op in MATCH_OPS
    )
    if matched_length == 0:
        raise ValueError(
            f"cigar at line {line_number} contains no M, =, or X operations: {cigar!r}"
        )

    return matched_length


def process_insertion_positions(input_file, output_file):
    """Read the input TSV, calculate insertion_pos, and write the results."""
    required_columns = [
        "sequence_id",
        "target_chrom",
        "match_start_pos",
        "cigar",
        "aligned_strand",
    ]

    with open(input_file, "r", encoding="utf-8", newline="") as input_handle, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input file is missing a header")

        missing_columns = [
            column for column in required_columns if column not in reader.fieldnames
        ]
        if missing_columns:
            raise ValueError(
                "Missing required fields: " + ", ".join(missing_columns)
            )

        writer = csv.writer(output_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("sequence_id", "target_chrom", "insertion_pos", "strand"))

        for line_number, row in enumerate(reader, start=2):
            try:
                match_start_pos = int(row["match_start_pos"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"match_start_pos is not an integer at line {line_number}: "
                    f"{row['match_start_pos']!r}"
                ) from error

            if match_start_pos < 1:
                raise ValueError(
                    f"match_start_pos must be greater than or equal to 1 at line {line_number}"
                )

            matched_len = matched_length_from_cigar(
                row["cigar"].strip(),
                line_number,
            )
            strand = row["aligned_strand"].strip().lower()

            if strand == "positive":
                insertion_pos = match_start_pos
            elif strand == "negative":
                insertion_pos = match_start_pos + matched_len - 1
            else:
                raise ValueError(
                    f"Unknown aligned_strand value at line {line_number}: "
                    f"{row['aligned_strand']!r}"
                )

            writer.writerow(
                (
                    row["sequence_id"],
                    row["target_chrom"],
                    insertion_pos,
                    strand,
                )
            )


def main():
    """Check command-line arguments and calculate insertion positions."""
    if len(sys.argv) != 3:
        program = sys.argv[0] if sys.argv else "get_insertion_pos.py"
        print(f"Usage: python {program} <input.tsv> <output.tsv>", file=sys.stderr)
        return 2

    try:
        process_insertion_positions(sys.argv[1], sys.argv[2])
    except (OSError, ValueError, csv.Error) as error:
        print(f"Processing failed: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())