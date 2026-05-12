#!/usr/bin/env python3
"""
Filter mismatched adapter assignments from paired trimmed FASTQ files.

Reads both input files entirely into memory in a single read() call to
minimize network filesystem round trips, then processes records in memory.
Pairs where R1 and R2 have matching adapter assignments are written to
.matched.1.fastq / .matched.2.fastq. Discarded pairs are written to
.discarded.1.fastq / .discarded.2.fastq. A summary CSV is written with
counts of concordant and discarded pairs.

Usage:
    filter_mismatched_adapters.py -r1 TRIMMED_R1 -r2 TRIMMED_R2 -s SAMPLE -o OUT_DIR
"""

import re
import sys
import argparse
import csv
from pathlib import Path

# Compile regex once at module level
ADAPTER_RE = re.compile(r'adapter=([^=\s]+)=')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter paired FASTQ records with mismatched adapter assignments."
    )
    parser.add_argument("-r1", "--read1", required=True,
                        help="Trimmed R1 FASTQ file.")
    parser.add_argument("-r2", "--read2", required=True,
                        help="Trimmed R2 FASTQ file.")
    parser.add_argument("-s", "--sample", required=True,
                        help="Sample name (used for output filenames and log).")
    parser.add_argument("-o", "--out_dir", required=True,
                        help="Output directory for matched and discarded FASTQ files.")
    return parser.parse_args()


def extract_adapter(header):
    """Extract adapter name from header field adapter=<name>=<sample>."""
    m = ADAPTER_RE.search(header)
    return m.group(1) if m else None


def load_records(filepath):
    """
    Read entire file in one call and split into a list of
    (header, seq, plus, qual) tuples.
    """
    with open(filepath, 'r', buffering=8*1024*1024) as fh:
        content = fh.read()

    lines = content.splitlines()
    if len(lines) % 4 != 0:
        print(f"WARNING: {filepath} has {len(lines)} lines, not a multiple of 4.",
              file=sys.stderr)

    records = []
    for i in range(0, len(lines) - 3, 4):
        records.append((lines[i], lines[i+1], lines[i+2], lines[i+3]))
    return records


def write_records(records, filepath):
    """Write a list of records to a file in a single write() call."""
    with open(filepath, 'w', buffering=8*1024*1024) as fh:
        fh.write('\n'.join(
            line
            for rec in records
            for line in rec
        ) + ('\n' if records else ''))


def main():
    args = parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    matched_r1  = out_dir / (args.sample + ".matched.1.fastq")
    matched_r2  = out_dir / (args.sample + ".matched.2.fastq")
    discard_r1  = out_dir / (args.sample + ".discarded.1.fastq")
    discard_r2  = out_dir / (args.sample + ".discarded.2.fastq")
    summary_csv = out_dir / (args.sample + ".adapter_filter.csv")

    # Load both files entirely into memory in one read() call each
    print(f"Loading {args.read1} ...")
    r1_records = load_records(args.read1)
    print(f"Loading {args.read2} ...")
    r2_records = load_records(args.read2)

    if len(r1_records) != len(r2_records):
        print(
            f"WARNING: R1 has {len(r1_records):,} records but "
            f"R2 has {len(r2_records):,} records. "
            f"Processing only paired records.",
            file=sys.stderr
        )

    total      = min(len(r1_records), len(r2_records))
    matched_r1_out = []
    matched_r2_out = []
    discard_r1_out = []
    discard_r2_out = []

    for rec1, rec2 in zip(r1_records, r2_records):
        adapter1 = extract_adapter(rec1[0])
        adapter2 = extract_adapter(rec2[0])

        if adapter1 is None:
            print(f"WARNING: Could not parse adapter from R1 header: {rec1[0]}",
                  file=sys.stderr)
        if adapter2 is None:
            print(f"WARNING: Could not parse adapter from R2 header: {rec2[0]}",
                  file=sys.stderr)

        if adapter1 is not None and adapter2 is not None and adapter1 == adapter2:
            matched_r1_out.append(rec1)
            matched_r2_out.append(rec2)
        else:
            discard_r1_out.append(rec1)
            discard_r2_out.append(rec2)

    concordant = len(matched_r1_out)
    discarded  = len(discard_r1_out)
    discard_pct = round(100.0 * discarded / total, 2) if total > 0 else 0.0

    # Write all output files in single write() calls
    print("Writing matched records ...")
    write_records(matched_r1_out, matched_r1)
    write_records(matched_r2_out, matched_r2)

    print("Writing discarded records ...")
    write_records(discard_r1_out, discard_r1)
    write_records(discard_r2_out, discard_r2)

    with open(summary_csv, 'w', newline='') as csv_fh:
        writer = csv.writer(csv_fh)
        writer.writerow(["sample", "total_pairs", "concordant_pairs",
                         "discarded_pairs", "discarded_pct"])
        writer.writerow([args.sample, total, concordant, discarded, discard_pct])

    print(f"{args.sample}: {total:,} total pairs, "
          f"{concordant:,} concordant, "
          f"{discarded:,} discarded ({discard_pct}%)")


if __name__ == "__main__":
    main()
