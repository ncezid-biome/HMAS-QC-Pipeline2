#!/bin/bash

cd ..
nextflow run hmas2.nf -profile test,git_action,singularity
echo "Looking for output in: $(pwd)"
find . -maxdepth 1 -type d | sort
NXF_EXIT=$?

if [ $NXF_EXIT -ne 0 ]; then
  echo "Nextflow pipeline failed"
  exit 1
fi

echo "Looking for output in: $(pwd)"
find . -maxdepth 1 -type d | sort
latest_testoutput=$(find . -maxdepth 1 -type d -name 'test_output*' | sort -r | head -n 1)

if [ -z "$latest_testoutput" ]; then
  echo "ERROR: Could not find test output directory"
  exit 1
fi

generated_csv="$latest_testoutput/report_sorted.csv"
expected_csv="test_data/report_ref_sorted.csv"

echo "Looking for output in: $(pwd)"
find . -maxdepth 1 -type d | sort
report_file=$(find "$latest_testoutput" -type f -name 'report*.csv')

if [ -z "$report_file" ]; then
  echo "ERROR: Could not find report CSV in $latest_testoutput"
  exit 1
fi

sort "$report_file" > "$generated_csv"
sort test_data/report_ref.csv > "$expected_csv"

if ! diff -q "$generated_csv" "$expected_csv" > /dev/null; then
  echo "WARNING ! *** CSV files differ ***"
  diff "$generated_csv" "$expected_csv"
else
  echo "PASSED ! CSV files match"
fi

rm "$generated_csv" "$expected_csv"