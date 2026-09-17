#!/usr/bin/env python3
import argparse
import re


FILTER_LENGTH = 50


def read_seq_len_from_fasta(fasta_file):
    with open(fasta_file, "r", encoding="utf-8") as fasta_handle:
        for line_number, line in enumerate(fasta_handle, start=1):
            if line_number == 2:
                sequence = line.strip()
                return len(sequence)

    raise ValueError(f"Cannot read sequence from second line of {fasta_file}")


def revcomp(seq):
    comp = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }
    return "".join(comp.get(base, "N") for base in reversed(seq))


def filter_sam(sam_file, output_file, seq_len, filter_length=FILTER_LENGTH):
    if filter_length < 0:
        raise ValueError("filter_length cannot be negative")

    pattern_left_softclip = re.compile(rf"^[0-9]+S{seq_len}M$")
    pattern_right_softclip = re.compile(rf"^{seq_len}M[0-9]+S$")

    with open(sam_file, "r", encoding="utf-8") as input_handle, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as output_handle:
        for line in input_handle:
            line = line.rstrip("\n")
            if not line or line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) < 11:
                continue

            read_id = fields[0]
            cigar = fields[5]
            sequence = fields[9]

            if pattern_left_softclip.fullmatch(cigar):
                reverse_sequence = revcomp(sequence)
                matched_seq = reverse_sequence[:seq_len]
                unmatched_seq = reverse_sequence[seq_len:]
            elif pattern_right_softclip.fullmatch(cigar):
                matched_seq = sequence[:seq_len]
                unmatched_seq = sequence[seq_len:]
            else:
                continue

            if len(unmatched_seq) < filter_length:
                continue

            output_handle.write(
                f"{read_id}\t{cigar}\t{matched_seq}\t{unmatched_seq}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter deduped SAM by CIGAR and extract matched/unmatched sequences. "
            f"Minimum unmatched sequence length: {FILTER_LENGTH}."
        )
    )
    parser.add_argument(
        "-f",
        "--guide_fasta",
        required=True,
    )
    parser.add_argument(
        "-i",
        "--input_sam",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_tsv",
        required=True,
    )
    args = parser.parse_args()

    seq_len = read_seq_len_from_fasta(args.guide_fasta)
    filter_sam(
        sam_file=args.input_sam,
        output_file=args.output_tsv,
        seq_len=seq_len,
    )


if __name__ == "__main__":
    main()
