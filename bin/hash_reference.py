#!/usr/bin/env python3
"""
hash_reference.py

Reads a multi-FASTA reference file, replaces each header line with the MD5
hash of that header, and writes:
  1. A new FASTA file with hashed headers (sequence unchanged, wrapped at 60 chars)
  2. A translation table (TSV) mapping original_id -> md5_hash

Hash collisions (extremely unlikely) are disambiguated with a numeric suffix.

Usage:
    hash_reference.py --input reference.fasta \\
                       --out-fasta reference.hashed.fasta \\
                       --out-table id_to_md5.tsv
"""

import argparse
import hashlib
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Hash FASTA headers and emit a translation table.")
    parser.add_argument("--input", required=True, help="Input multi-FASTA reference file")
    parser.add_argument("--out-fasta", required=True, help="Output FASTA file with MD5-hashed headers")
    parser.add_argument("--out-table", required=True, help="Output TSV translation table (original_id, md5_hash)")
    parser.add_argument("--line-width", type=int, default=60, help="Sequence line wrap width (default: 60)")
    return parser.parse_args()


def md5_of(text):
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def read_fasta_records(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_chunks = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if header is not None:
            yield header, "".join(seq_chunks)


def main():
    args = parse_args()

    records = list(read_fasta_records(args.input))
    if not records:
        sys.exit(f"ERROR: No FASTA records found in {args.input}")

    seen_hashes = {}  # candidate_hash -> original_header (for collision detection)

    with open(args.out_fasta, "w") as out_fa, open(args.out_table, "w") as out_tbl:
        out_tbl.write("original_id\tmd5_hash\n")

        for header, seq in records:
            # original_id = first whitespace-delimited token of the header
            original_id = header.split()[0] if header else ""
            h = md5_of(header)

            # Disambiguate on the rare chance of a collision
            suffix = 0
            candidate = h
            while candidate in seen_hashes and seen_hashes[candidate] != header:
                suffix += 1
                candidate = f"{h}_{suffix}"
            seen_hashes[candidate] = header

            out_fa.write(f">{candidate}\n")
            for i in range(0, len(seq), args.line_width):
                out_fa.write(seq[i:i + args.line_width] + "\n")

            out_tbl.write(f"{original_id}\t{candidate}\n")

    print(f"Processed {len(records)} reference sequences.")


if __name__ == "__main__":
    main()
