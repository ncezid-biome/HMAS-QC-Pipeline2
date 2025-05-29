#!/usr/bin/env python

import argparse
import shlex

def parse_kv_pairs(input_string):
    tokens = shlex.split(input_string)
    pairs = []
    i = 0
    while i < len(tokens) - 1:
        if tokens[i].startswith('--'):
            key = tokens[i]
            value = tokens[i + 1]
            if not value.startswith('--'):  # assume it's a value
                pairs.append((key, value))
                i += 2
                continue
        i += 1
    return pairs

def main():
    parser = argparse.ArgumentParser(description='Generate HTML block for MultiQC custom content.')
    parser.add_argument('--cmdline', required=True, help='The full Nextflow workflow.commandLine string')
    parser.add_argument('--params_str', required=False, help='Other Nextflow param key-value string')
    parser.add_argument('--output', required=True, help='Output file')

    args = parser.parse_args()

    # Parse both commandLine and params_str into key/value pairs
    cmdline_pairs = parse_kv_pairs(args.cmdline)
    param_pairs = parse_kv_pairs(args.params_str) if args.params_str else []

    # all_pairs = cmdline_pairs + param_pairs
    all_pairs = param_pairs + cmdline_pairs

    # remove duplicates by keeping the last value for each key
    unique_pairs = {}
    for k, v in all_pairs:
        unique_pairs[k] = v

    with open(args.output, 'w') as f:
        f.write("# id: 'cli-html'\n")
        f.write("# section_name: 'Pipeline Arguments'\n")
        f.write("# description: 'This section lists some of the main arguments used'\n")
        f.write("# plot_type: 'html'\n\n")
        f.write("<dl class=dl-horizontal>\n")
        # for key, value in all_pairs:
        for key, value in unique_pairs.items():
            f.write(f"  <dt>{key}</dt>\n")
            f.write(f"  <dd>{value}</dd>\n")
        f.write("</dl>\n")

if __name__ == "__main__":
    main()
