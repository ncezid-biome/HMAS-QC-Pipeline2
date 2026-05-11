#!/usr/bin/env python3
"""
Bin FASTQ records by adapter name found in the read header.

Pass 1: Scan headers only to build a dict mapping read_id -> adapter,
        independently for R1 and R2. Performs concordance checks.
Pass 2: For each adapter, stream through the input file once and write
        only matching records to the corresponding output file.
        Only 2 file handles are open at any time (input + current output).

Usage:
    python bin_fastq_by_adapter.py -i INPUT_DIR -o OUTPUT_DIR
"""

import os
import re
import sys
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Bin paired FASTQ records by adapter= field in read headers."
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory containing paired trimmed FASTQ files.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory where output subdirectories and files will be created.")
    return parser.parse_args()


def extract_adapter(header):
    m = re.search(r'adapter=([^=\s]+)=', header)
    return m.group(1) if m else None


def extract_read_id(header):
    return header.lstrip('@').split()[0]


def pass1_build_index(filepath):
    read_index = {}
    adapters = set()
    with open(filepath, 'r') as fh:
        while True:
            header = fh.readline().rstrip('\n')
            if not header:
                break
            if not header.startswith('@'):
                print(f"WARNING: Unexpected header line: {header}", file=sys.stderr)
                for _ in range(3):
                    fh.readline()
                continue
            for _ in range(3):
                fh.readline()
            read_id = extract_read_id(header)
            adapter = extract_adapter(header)
            if adapter is None:
                print(f"WARNING: Could not parse adapter from header: {header}", file=sys.stderr)
                adapter = "unknown"
            read_index[read_id] = adapter
            adapters.add(adapter)
    return read_index, adapters


def check_concordance(r1_index, r1_adapters, r2_index, r2_adapters):
    adapters = r1_adapters | r2_adapters

    if r1_adapters == r2_adapters:
        print("  Adapter sets match between R1 and R2.")
    else:
        only_r1 = r1_adapters - r2_adapters
        only_r2 = r2_adapters - r1_adapters
        if only_r1:
            print(f"  WARNING: Adapters in R1 only: {only_r1}", file=sys.stderr)
        if only_r2:
            print(f"  WARNING: Adapters in R2 only: {only_r2}", file=sys.stderr)

    r1_ids = set(r1_index.keys())
    r2_ids = set(r2_index.keys())
    if r1_ids == r2_ids:
        print("  Read IDs match between R1 and R2.")
    else:
        only_r1_ids = r1_ids - r2_ids
        only_r2_ids = r2_ids - r1_ids
        if only_r1_ids:
            print(f"  WARNING: {len(only_r1_ids):,} read ID(s) in R1 not found in R2. "
                  f"First 5: {list(only_r1_ids)[:5]}", file=sys.stderr)
        if only_r2_ids:
            print(f"  WARNING: {len(only_r2_ids):,} read ID(s) in R2 not found in R1. "
                  f"First 5: {list(only_r2_ids)[:5]}", file=sys.stderr)

    shared_ids = r1_ids & r2_ids
    mismatched = {rid for rid in shared_ids if r1_index[rid] != r2_index[rid]}
    if mismatched:
        print(f"  WARNING: {len(mismatched):,} read ID(s) have different adapter "
              f"assignments in R1 vs R2. First 5:", file=sys.stderr)
        for rid in list(mismatched)[:5]:
            print(f"    {rid}: R1={r1_index[rid]}, R2={r2_index[rid]}", file=sys.stderr)
    else:
        print("  Adapter assignments are consistent for all shared read IDs.")

    return adapters


def pass2_write_by_adapter(filepath, read_index, adapters, out_root, base, read_num):
    counts = {adapter: 0 for adapter in adapters}
    for adapter in sorted(adapters):
        out_path = out_root / (base + "." + str(read_num) + "." + adapter + ".fastq")
        with open(filepath, 'r') as in_fh, open(out_path, 'w') as out_fh:
            while True:
                header = in_fh.readline().rstrip('\n')
                seq    = in_fh.readline().rstrip('\n')
                plus   = in_fh.readline().rstrip('\n')
                qual   = in_fh.readline().rstrip('\n')
                if not header:
                    break
                if not header.startswith('@'):
                    print(f"WARNING: Unexpected header line: {header}", file=sys.stderr)
                    continue
                read_id = extract_read_id(header)
                record_adapter = read_index.get(read_id, "unknown")
                if record_adapter == adapter:
                    out_fh.write(header + "\n" + seq + "\n" + plus + "\n" + qual + "\n")
                    counts[adapter] += 1
    return counts


def find_pairs(input_dir):
    input_path = Path(input_dir)
    r1_files = sorted(input_path.rglob("*.trimmed.1.fastq"))
    if not r1_files:
        print(f"ERROR: No *.trimmed.1.fastq files found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    pairs = []
    for r1 in r1_files:
        stem = r1.name
        common_root = stem[:-len(".1.fastq")]
        r2 = r1.parent / (common_root + ".2.fastq")
        if not r2.exists():
            print(f"WARNING: No R2 partner found for {r1.name}, skipping.", file=sys.stderr)
            continue
        pairs.append((r1, r2, common_root))
    return pairs


def process_pair(r1_path, r2_path, common_root, output_dir):
    out_root = Path(output_dir) / common_root
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"  Pass 1: indexing headers from {r1_path.name} ...")
    r1_index, r1_adapters = pass1_build_index(r1_path)
    print(f"  Found {len(r1_index):,} records across {len(r1_adapters)} adapter(s).")

    print(f"  Pass 1: indexing headers from {r2_path.name} ...")
    r2_index, r2_adapters = pass1_build_index(r2_path)
    print(f"  Found {len(r2_index):,} records across {len(r2_adapters)} adapter(s).")

    adapters = check_concordance(r1_index, r1_adapters, r2_index, r2_adapters)

    base = r1_path.name[:-len(".1.fastq")]

    print(f"  Pass 2: writing R1 records ({len(adapters)} adapters) ...")
    r1_counts = pass2_write_by_adapter(r1_path, r1_index, adapters, out_root, base, 1)

    print(f"  Pass 2: writing R2 records ({len(adapters)} adapters) ...")
    r2_counts = pass2_write_by_adapter(r2_path, r2_index, adapters, out_root, base, 2)

    print("  Results:")
    for adapter in sorted(adapters):
        print(f"    [{adapter}] R1: {r1_counts[adapter]:,}, R2: {r2_counts[adapter]:,} records")


def main():
    args = parse_args()
    if not os.path.isdir(args.input_dir):
        print(f"ERROR: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)
    pairs = find_pairs(args.input_dir)
    print(f"Found {len(pairs)} FASTQ pair(s) to process.\n")
    for r1, r2, common_root in pairs:
        print(f"Processing: {common_root}")
        process_pair(r1, r2, common_root, args.output_dir)
        print()
    print("Done.")


if __name__ == "__main__":
    main()
