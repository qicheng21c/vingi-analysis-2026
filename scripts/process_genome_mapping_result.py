#!/usr/bin/env python3
import re
import sys


READ_CONSUMING_OPS = {"M", "I", "S", "=", "X"}
REFERENCE_CONSUMING_OPS = {"M", "D", "N", "=", "X"}
CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def reverse_complement(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1].upper()


def parse_cigar(cigar):
    """Fully parse the CIGAR string and return None if invalid."""
    ops = [(int(length), op) for length, op in CIGAR_PATTERN.findall(cigar)]

    if not ops or any(length <= 0 for length, _ in ops):
        return None

    reconstructed = "".join(f"{length}{op}" for length, op in ops)
    if reconstructed != cigar:
        return None

    return ops


def read_consuming_length(ops):
    """Calculate the total read-consuming length in the CIGAR string."""
    return sum(length for length, op in ops if op in READ_CONSUMING_OPS)


def reference_consuming_length(ops):
    """Calculate the total reference-consuming length in the CIGAR string."""
    return sum(length for length, op in ops if op in REFERENCE_CONSUMING_OPS)


def five_prime_softclip_length(ops, flag):
    """Return the length of the soft-clipped region at the 5' end of the original read."""
    oriented_ops = ops if flag == 0 else reversed(ops)
    operation_iterator = iter(oriented_ops)

    for length, op in operation_iterator:
        if op == "H":
            continue
        if op != "S":
            return 0

        softclip_length = length
        for next_length, next_op in operation_iterator:
            if next_op != "S":
                break
            softclip_length += next_length
        return softclip_length

    return 0


def split_oriented_read(seq, softclip_length):
    junction_seq = seq[:softclip_length]
    genome_seq = seq[softclip_length:]
    return junction_seq, genome_seq


def write_record(output_handle, fields):
    output_handle.write("\t".join(str(field) for field in fields) + "\n")


def process_sam(input_sam, output_tsv):
    header = (
        "sequence_id",
        "target_chrom",
        "match_start_pos",
        "cigar",
        "junction_seq",
        "genome_seq",
        "aligned_strand",
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
                pos = int(fields[3])
            except ValueError:
                continue

            if flag not in {0, 16}:
                continue

            sequence_id = fields[0]
            target_chrom = fields[2]
            cigar = fields[5]
            seq = fields[9]

            if cigar == "*" or seq == "*":
                continue

            ops = parse_cigar(cigar)
            if ops is None:
                continue

            if read_consuming_length(ops) != len(seq):
                continue

            reference_length = reference_consuming_length(ops)
            if reference_length == 0:
                continue

            softclip_length = five_prime_softclip_length(ops, flag)

            if flag == 0:
                aligned_strand = "positive"
                match_start_pos = pos
                oriented_seq = seq
            else:
                aligned_strand = "negative"
                match_start_pos = pos + reference_length - 1
                oriented_seq = reverse_complement(seq)

            junction_seq, genome_seq = split_oriented_read(
                oriented_seq,
                softclip_length,
            )

            write_record(
                output_handle,
                (
                    sequence_id,
                    target_chrom,
                    match_start_pos,
                    cigar,
                    junction_seq,
                    genome_seq,
                    aligned_strand,
                ),
            )


def main():
    if len(sys.argv) != 3:
        program = sys.argv[0] if sys.argv else "process_genome_mapping_result.py"
        print(f"Usage: {program} <input.sam> <output.tsv>", file=sys.stderr)
        return 2

    try:
        process_sam(sys.argv[1], sys.argv[2])
    except OSError as error:
        print(f"Failed: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
