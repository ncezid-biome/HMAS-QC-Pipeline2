#!/usr/bin/env python3
"""
Bin FASTQ records by adapter name found in the read header.

Pass 1: Scan headers only to build a dict mapping read_id -> adapter,
        independently for R1 and R2. Performs concordance checks.
Pass 2: Read each file once into memory, bin records by adapter in memory,
        then write each adapter's records sequentially with only 2 file
        handles open at a time (never more than one output file at once).

Usage:
    python bin_fastq_by_adapter.py -r1 R1.fastq -r2 R2.fastq -o OUTPUT_DIR
"""

import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict


# Compile once at module level
ADAPTER_RE = re.compile(r'adapter=([^=\s]+)=')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Bin paired FASTQ records by adapter= field in read headers."
    )
    parser.add_argument("-r1", "--read1", required=True,
                        help="R1 matched FASTQ file.")
    parser.add_argument("-r2", "--read2", required=True,
                        help="R2 matched FASTQ file.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory where output files will be created.")
    return parser.parse_args()


def extract_adapter(header):
    m = ADAPTER_RE.search(header)
    return m.group(1) if m else None


def extract_read_id(header):
    return header.lstrip('@').split()[0]


def pass1_build_index(filepath):
    """
    First pass: scan headers only, building a dict of read_id -> adapter.
    Reads entire file at once to minimize network filesystem round trips.
    Returns (read_index, adapters).
    """
    read_index = {}
    adapters = set()

    with open(filepath, 'r', buffering=8*1024*1024) as fh:
        content = fh.read()

    lines = content.splitlines()
    for i in range(0, len(lines), 4):
        if i >= len(lines) or not lines[i].startswith('@'):
            continue
        header = lines[i]
        read_id = extract_read_id(header)
        adapter = extract_adapter(header)
        if adapter is None:
            print(f"WARNING: Could not parse adapter from header: {header}",
                  file=sys.stderr)
            adapter = "unknown"
        read_index[read_id] = adapter
        adapters.add(adapter)

    return read_index, adapters


def check_concordance(r1_index, r1_adapters, r2_index, r2_adapters):
    """
    Check R1 and R2 for concordance of adapter sets, read IDs,
    and adapter assignments for shared read IDs.
    Returns the union of adapter names across both files.
    """
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
            print(f"    {rid}: R1={r1_index[rid]}, R2={r2_index[rid]}",
                  file=sys.stderr)
    else:
        print("  Adapter assignments are consistent for all shared read IDs.")

    return adapters


def pass2_bin_and_write(filepath, read_index, adapters, out_dir, base, read_num):
    """
    Second pass: read entire file into memory once, bin records by adapter,
    then write each adapter's records to its output file sequentially.
    Only one output file handle is open at a time.
    Returns a dict of {adapter: record_count}.
    """
    # Read entire file in one call
    with open(filepath, 'r', buffering=8*1024*1024) as fh:
        content = fh.read()

    lines = content.splitlines()

    # Bin records in memory: {adapter: [list of joined record strings]}
    bins = defaultdict(list)
    for i in range(0, len(lines) - 3, 4):
        if not lines[i].startswith('@'):
            print(f"WARNING: Unexpected header line: {lines[i]}", file=sys.stderr)
            continue
        header = lines[i]
        record_str = header + "\n" + lines[i+1] + "\n" + lines[i+2] + "\n" + lines[i+3]
        read_id = extract_read_id(header)
        adapter = read_index.get(read_id, "unknown")
        bins[adapter].append(record_str)

    # Write each adapter's records sequentially, one file at a time
    counts = {}
    for adapter in sorted(adapters):
        out_path = out_dir / (base + "." + str(read_num) + "." + adapter + ".fastq")
        records = bins.get(adapter, [])
        with open(out_path, 'w', buffering=8*1024*1024) as fh:
            if records:
                fh.write('\n'.join(records) + '\n')
        counts[adapter] = len(records)

    return counts


def main():
    args = parse_args()

    r1_path = Path(args.read1)
    r2_path = Path(args.read2)
    out_dir = Path(args.output_dir)

    if not r1_path.exists():
        print(f"ERROR: R1 file not found: {r1_path}", file=sys.stderr)
        sys.exit(1)
    if not r2_path.exists():
        print(f"ERROR: R2 file not found: {r2_path}", file=sys.stderr)
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    # Derive base name from R1 filename by stripping .1.fastq suffix
    base = r1_path.name[:-len(".1.fastq")]

    print(f"Processing: {base}")

    # Pass 1: build independent indices for R1 and R2
    print(f"  Pass 1: indexing headers from {r1_path.name} ...")
    r1_index, r1_adapters = pass1_build_index(r1_path)
    print(f"  Found {len(r1_index):,} records across {len(r1_adapters)} adapter(s).")

    print(f"  Pass 1: indexing headers from {r2_path.name} ...")
    r2_index, r2_adapters = pass1_build_index(r2_path)
    print(f"  Found {len(r2_index):,} records across {len(r2_adapters)} adapter(s).")

    adapters = check_concordance(r1_index, r1_adapters, r2_index, r2_adapters)

    # Pass 2: bin into memory and write sequentially
    print(f"  Pass 2: binning and writing R1 records ({len(adapters)} adapters) ...")
    r1_counts = pass2_bin_and_write(r1_path, r1_index, adapters, out_dir, base, 1)

    print(f"  Pass 2: binning and writing R2 records ({len(adapters)} adapters) ...")
    r2_counts = pass2_bin_and_write(r2_path, r2_index, adapters, out_dir, base, 2)

    print("  Results:")
    for adapter in sorted(adapters):
        print(f"    [{adapter}] R1: {r1_counts[adapter]:,}, R2: {r2_counts[adapter]:,} records")

    print("Done.")


if __name__ == "__main__":
    main()
