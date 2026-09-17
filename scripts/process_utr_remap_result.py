#!/usr/bin/env python3
import re
import sys


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")
READ_CONSUMING_OPS = {"M", "I", "S", "=", "X"}
SCAR_START_OPS = {"I", "D", "S"}


def parse_cigar(cigar):
    """Fully parse the CIGAR string and return None if the format is invalid."""
    ops = [(int(length), op) for length, op in CIGAR_PATTERN.findall(cigar)]

    if not ops or any(length <= 0 for length, _ in ops):
        return None

    reconstructed = "".join(f"{length}{op}" for length, op in ops)
    if reconstructed != cigar:
        return None

    return ops


def find_scar_boundary(ops):
    """
    Return the operation index of the first I, D, or S.

    Return the length of ops if the CIGAR contains only M operations,
    or None if the CIGAR structure is invalid.
    """
    if not ops or ops[0][1] != "M":
        return None

    for index, (_, op) in enumerate(ops):
        if op in SCAR_START_OPS:
            # Only successfully aligned M operations are allowed before the split boundary.
            if all(previous_op == "M" for _, previous_op in ops[:index]):
                return index
            return None

        if op != "M":
            return None

    # If the CIGAR contains only M operations, the split boundary is at the end
    # of the CIGAR and scar_seq is empty.
    return len(ops)


def read_consuming_length(ops):
    """Calculate the total read-consuming length in the CIGAR string."""
    return sum(length for length, op in ops if op in READ_CONSUMING_OPS)


def split_seq_by_cigar(seq, ops, boundary):
    """Split the sequence into left-side matched_seq and right-side scar_seq at the CIGAR boundary."""
    seq_pos = 0
    matched_parts = []
    scar_parts = []

    for index, (length, op) in enumerate(ops):
        if op not in READ_CONSUMING_OPS:
            continue

        sequence_part = seq[seq_pos:seq_pos + length]
        if index < boundary:
            matched_parts.append(sequence_part)
        else:
            scar_parts.append(sequence_part)
        seq_pos += length

    return "".join(matched_parts), "".join(scar_parts)


def write_record(output_handle, fields):
    """Write one TSV record to the output file."""
    output_handle.write("\t".join(str(field) for field in fields) + "\n")


def process_sam(input_sam, output_tsv):
    """Read the SAM file, split valid records, and write them to a TSV file."""
    header = (
        "sequence_id",
        "target_seq",
        "match_start_pos",
        "cigar",
        "matched_seq",
        "scar_seq",
    )

    with open(input_sam, "r", encoding="utf-8") as input_handle, open(
        output_tsv, "w", encoding="utf-8", newline=""
    ) as output_handle:
        write_record(output_handle, header)

        for line in input_handle:
            line = line.rstrip("\r\n")
            if not line or line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) < 11:
                continue

            try:
                flag = int(fields[1])
                match_start_pos = int(fields[3])
            except ValueError:
                continue

            if flag not in {0, 16}:
                continue

            sequence_id = fields[0]
            target_seq = fields[2]
            cigar = fields[5]
            seq = fields[9]

            if cigar == "*" or seq == "*":
                continue

            ops = parse_cigar(cigar)
            boundary = find_scar_boundary(ops)
            if boundary is None:
                continue

            if read_consuming_length(ops) != len(seq):
                continue

            matched_seq, scar_seq = split_seq_by_cigar(seq, ops, boundary)

            write_record(
                output_handle,
                (
                    sequence_id,
                    target_seq,
                    match_start_pos,
                    cigar,
                    matched_seq,
                    scar_seq,
                ),
            )


def main():
    """Parse command-line arguments and run the conversion."""
    if len(sys.argv) != 3:
        program = sys.argv[0] if sys.argv else "process_utr_remap_result.py"
        print(f"Usage: {program} <input.sam> <output.tsv>", file=sys.stderr)
        return 2

    try:
        process_sam(sys.argv[1], sys.argv[2])
    except OSError as error:
        print(f"File processing failed: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())