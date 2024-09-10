#!/bin/bash

# Run Nextflow pipeline
cd ..
nextflow run hmas2.nf -profile test

latest_testoutput=$(ls -1t test_output* | head -n 1 | sed 's/[,;:]//g')

# Define file paths
generated_csv="$latest_testoutput/report_sorted.csv"
expected_csv="test_data/report_ref_sorted.csv"

# Check if the pipeline run was successful
if [ $? -ne 0 ]; then
  echo "Nextflow pipeline failed"
  exit 1
fi

report_file=$(find $latest_testoutput -type f -name 'report*.csv')

# Sort the generated CSV file and the expected CSV file
sort "$report_file" > "$generated_csv"
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
