#!/usr/bin/env python3
import csv
import re
import sys
from collections import defaultdict


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
REQUIRED_FIELDS = {
    "sequence_id",
    "target_seq",
    "match_start_pos",
    "cigar",
    "matched_seq",
    "scar_seq",
}


def parse_cigar(cigar, line_number):
    """Fully parse the CIGAR string and raise a ValueError with the input line number if the format is invalid."""
    ops = [(int(length), op) for length, op in CIGAR_PATTERN.findall(cigar)]

    if not ops or any(length <= 0 for length, _ in ops):
        raise ValueError(
            f"Invalid CIGAR at line {line_number}: {cigar!r}"
        )

    reconstructed = "".join(f"{length}{op}" for length, op in ops)
    if reconstructed != cigar:
        raise ValueError(
            f"Incomplete CIGAR format or unknown operation at line {line_number}: {cigar!r}"
        )

    return ops


def leading_match_length(ops):
    """Return the total length of consecutive M operations before the first non-M operation in the CIGAR."""
    match_length = 0

    for length, op in ops:
        if op != "M":
            break
        match_length += length

    return match_length


def calculate_map_depth(input_file):
    """Read the input TSV and return the difference counts and maximum covered position."""
    depth_changes = defaultdict(int)
    maximum_position = 0
    target_sequences = set()

    with open(input_file, "r", encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError("Input file is missing a header")

        missing_fields = REQUIRED_FIELDS.difference(reader.fieldnames)
        if missing_fields:
            missing_text = ", ".join(sorted(missing_fields))
            raise ValueError(f"Input file is missing required fields: {missing_text}")

        for line_number, row in enumerate(reader, start=2):
            target_seq = row["target_seq"].strip()
            if not target_seq:
                raise ValueError(f"target_seq is empty at line {line_number}")
            target_sequences.add(target_seq)

            try:
                match_start_pos = int(row["match_start_pos"])
            except ValueError as error:
                raise ValueError(
                    f"match_start_pos is not an integer at line {line_number}: "
                    f"{row['match_start_pos']!r}"
                ) from error

            if match_start_pos < 1:
                raise ValueError(
                    f"match_start_pos must be greater than or equal to 1 at line {line_number}"
                )

            ops = parse_cigar(row["cigar"].strip(), line_number)
            match_length = leading_match_length(ops)
            if match_length == 0:
                continue

            match_end_pos = match_start_pos + match_length - 1
            depth_changes[match_start_pos] += 1
            depth_changes[match_end_pos + 1] -= 1
            maximum_position = max(maximum_position, match_end_pos)

    if len(target_sequences) > 1:
        targets = ", ".join(sorted(target_sequences))
        raise ValueError(
            "Input file contains multiple target_seq values and cannot be merged "
            "in output without distinguishing references: "
            f"{targets}"
        )

    return depth_changes, maximum_position


def write_map_depth(output_file, depth_changes, maximum_position):
    """Write the map depth for each reference position based on the difference counts."""
    current_depth = 0

    with open(output_file, "w", encoding="utf-8", newline="") as output_handle:
        output_handle.write("position\tmap_depth\n")

        for position in range(1, maximum_position + 1):
            current_depth += depth_changes.get(position, 0)
            output_handle.write(f"{position}\t{current_depth}\n")


def process_map_depth(input_file, output_file):
    """Calculate and write the UTR3 reference map depth."""
    depth_changes, maximum_position = calculate_map_depth(input_file)
    write_map_depth(output_file, depth_changes, maximum_position)


def main():
    """Parse command-line arguments and perform map depth calculation."""
    if len(sys.argv) != 3:
        program = sys.argv[0] if sys.argv else "utr3_map_depth_count.py"
        print(f"Usage: {program} <input> <output>", file=sys.stderr)
        return 2

    try:
        process_map_depth(sys.argv[1], sys.argv[2])
    except (OSError, ValueError) as error:
        print(f"Processing failed: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())