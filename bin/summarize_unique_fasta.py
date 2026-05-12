#!/usr/bin/env python3
"""
Summarize unique sequences per adapter from final.unique.fasta files.

Recursively searches an input directory for *.final.unique.fasta files.
For each file, produces a TSV alongside the fasta where each row contains:
    - adapter ID
    - total number of unique sequences for that adapter
    - abundance of each unique sequence in descending order

Sequence ID format expected:
    >READID=ADAPTER=SAMPLE;size=ABUNDANCE

Usage:
    python summarize_unique_fasta.py -i INPUT_DIR
"""

import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict


# Compile once at module level
SIZE_RE = re.compile(r';size=(\d+)')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize unique sequences per adapter from final.unique.fasta files."
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory to search recursively for *.final.unique.fasta files.")
    return parser.parse_args()


def extract_adapter_and_size(header):
    """
    Parse adapter name and size from a fasta header line.
    Expected format: >READID=ADAPTER=SAMPLE;size=ABUNDANCE
    Returns (adapter, size) or (None, None) if parsing fails.
    """
    # Strip leading '>'
    header = header.lstrip('>')

    # Extract size
    size_m = SIZE_RE.search(header)
    if size_m is None:
        return None, None
    size = int(size_m.group(1))

    # Extract adapter — second '='-delimited field
    fields = header.split('=')
    if len(fields) < 2:
        return None, None
    adapter = fields[1]

    return adapter, size


def process_fasta(fasta_path):
    """
    Read a fasta file entirely into memory and bin sequences by adapter.
    Returns a dict of {adapter: [size, size, ...]} with sizes in the
    order they appear in the file.
    """
    with open(fasta_path, 'r', buffering=8*1024*1024) as fh:
        content = fh.read()

    adapter_sizes = defaultdict(list)
    lines = content.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('>'):
            adapter, size = extract_adapter_and_size(line)
            if adapter is None:
                print(f"WARNING: Could not parse header: {line}", file=sys.stderr)
            else:
                adapter_sizes[adapter].append(size)
            i += 1
            # skip sequence line(s) until next header or EOF
            while i < len(lines) and not lines[i].startswith('>'):
                i += 1
        else:
            i += 1

    return adapter_sizes


def write_tsv(adapter_sizes, out_path):
    """
    Write TSV where each row is:
        adapter  total_unique  size1  size2  ... (descending order)
    """
    with open(out_path, 'w', buffering=8*1024*1024) as fh:
        # Header
        fh.write("adapter\ttotal_unique\tabundances\n")
        for adapter in sorted(adapter_sizes.keys()):
            sizes = sorted(adapter_sizes[adapter], reverse=True)
            total = len(sizes)
            sizes_str = '\t'.join(str(s) for s in sizes)
            fh.write(f"{adapter}\t{total}\t{sizes_str}\n")


def main():
    args = parse_args()
    input_path = Path(args.input_dir)

    if not input_path.is_dir():
        print(f"ERROR: Input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    fasta_files = sorted(input_path.rglob("*.final.unique.fasta"))
    if not fasta_files:
        print(f"ERROR: No *.final.unique.fasta files found in {args.input_dir}",
              file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(fasta_files)} fasta file(s) to process.")

    for fasta_path in fasta_files:
        print(f"Processing: {fasta_path.name} ...")
        adapter_sizes = process_fasta(fasta_path)

        out_path = fasta_path.parent / (fasta_path.stem + ".adapter_summary.tsv")
        write_tsv(adapter_sizes, out_path)

        total_adapters = len(adapter_sizes)
        total_seqs = sum(len(v) for v in adapter_sizes.values())
        print(f"  {total_adapters} adapters, {total_seqs:,} total unique sequences -> {out_path.name}")

    print("Done.")


if __name__ == "__main__":
    main()
