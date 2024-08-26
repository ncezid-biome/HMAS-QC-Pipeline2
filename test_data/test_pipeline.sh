#!/bin/bash

# Define file paths
generated_csv="test_output/report_sorted.csv"
expected_csv="test_data/report_ref_sorted.csv"

# Run Nextflow pipeline
cd ..
nextflow run hmas2.nf -profile test

# Check if the pipeline run was successful
if [ $? -ne 0 ]; then
  echo "Nextflow pipeline failed"
  exit 1
fi

# Sort the generated CSV file and the expected CSV file
sort test_output/report.csv > "$generated_csv"
sort test_data/report_ref.csv > "$expected_csv"

# Compare the sorted CSV files
if ! diff -q "$generated_csv" "$expected_csv" > /dev/null; then
  echo "WARNING ! *** CSV files differ ***"
  diff "$generated_csv" "$expected_csv"
else
  echo "PASSED ! CSV files match"
fi

# Clean up sorted temporary files
rm "$generated_csv" "$expected_csv"

# get rid of nasty ^H after running nextflow
stty erase ^H
