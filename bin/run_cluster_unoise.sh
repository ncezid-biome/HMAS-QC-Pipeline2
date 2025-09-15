#!/bin/bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <input_fasta> <oligos_file> <output_fasta> <minsize> <unoise_alpha>"
    exit 1
fi

input_fasta="$1"       # input fasta
oligos_file="$2"       # 4th col has OGxxxxprimerGroupX
output_file="$3"       # final output sorted fasta
minsize="$4"           # vsearch --minsize
unoise_alpha="$5"      # vsearch --unoise_alpha

if [[ ! -s "$input_fasta" ]]; then
    echo "⚠️ Input fasta is empty or missing. Exiting."
    exit 0
fi

# Temp directory (auto cleanup)
tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

concat_file="$tmpdir/concat.fasta"
> "$concat_file"   # truncate/init

# Loop through each pattern in column 4
awk '{print $4}' "$oligos_file" | tr -d '\r' | while read -r pattern; do

    subset="$tmpdir/${pattern}.fasta"
    awk -v pat="$pattern" '/^>/{keep=($0 ~ pat)} keep' "$input_fasta" > "$subset"

    if [[ -s "$subset" ]]; then
        out="$tmpdir/${pattern}.final.fasta"
        vsearch --cluster_unoise "$subset" \
                --minsize "$minsize" \
                --unoise_alpha "$unoise_alpha" \
                --centroids "$out" \
                --sizein --sizeout \
                --quiet
        cat "$out" >> "$concat_file"
    fi
done

# sort concatenated fasta by size= (descending) ===
if [[ -s "$concat_file" ]]; then
    vsearch --sortbysize "$concat_file" \
            --output   "$output_file" \
            --quiet
    echo "✅ All patterns processed. Concatenated and sorted fasta:  $output_file"
else
    echo "⚠️  No sequences were processed, nothing to concatenate or sort."
fi

