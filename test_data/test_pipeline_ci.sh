#!/bin/bash

cd ..
nextflow run hmas2.nf -profile test,git_action,singularity
NXF_EXIT=$?

if [ $NXF_EXIT -ne 0 ]; then
  echo "Nextflow pipeline failed"
  exit 1
fi

latest_testoutput=$(find . -maxdepth 1 -type d \( -name 'test_output*' -o -name 'hmas2_results*' \) | sort -r | head -n 1)

if [ -z "$latest_testoutput" ]; then
  echo "ERROR: Could not find test output directory"
  exit 1
fi

generated_csv="$latest_testoutput/report_sorted.csv"
expected_csv="test_data/report_ref_sorted.csv"
generated_tsv="$latest_testoutput/species_sorted.tsv"
expected_tsv="test_data/species_ref_sorted.tsv"

report_file=$(find "$latest_testoutput" -type f -name 'report*.csv')
species_file=$(find "$latest_testoutput" -type f -name 'species*.tsv')

if [ -z "$report_file" ]; then
  echo "ERROR: Could not find report CSV in $latest_testoutput"
  exit 1
fi

if [ -z "$species_file" ]; then
  echo "ERROR: Could not find species TSV in $latest_testoutput"
  exit 1
fi

sort "$report_file" > "$generated_csv"
sort test_data/report_ref.csv > "$expected_csv"
sort "$species_file" > "$generated_tsv"
sort test_data/species_ref.tsv > "$expected_tsv"

if ! diff -q "$generated_csv" "$expected_csv" > /dev/null; then
  echo "WARNING ! *** CSV files differ ***"
  diff "$generated_csv" "$expected_csv"
else
  echo "PASSED ! CSV files match"
fi

if ! diff -q "$generated_tsv" "$expected_tsv" > /dev/null; then
  echo "WARNING ! *** TSV files differ ***"
  diff "$generated_tsv" "$expected_tsv"
else
  echo "PASSED ! TSV files match"
fi

rm "$generated_csv" "$expected_csv" "$generated_tsv" "$expected_tsv"
