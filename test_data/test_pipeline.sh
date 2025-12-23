#!/usr/bin/env bash
set -euo pipefail

############################################
# Configuration (override via env if needed)
############################################

NF_PROFILE="${NF_PROFILE:-test}"
NF_PIPELINE="hmas2.nf"
TEST_DATA_DIR="test_data"
REF_CSV="$TEST_DATA_DIR/report_ref.csv"

############################################
# Helper functions
############################################

die() {
  echo "ERROR: $*" >&2
  exit 1
}

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

############################################
# Pre-flight checks
############################################

command_exists nextflow || die "Nextflow not found in PATH"

# Detect container runtime
if command_exists singularity; then
  CONTAINER_PROFILE="singularity"
elif command_exists apptainer; then
  CONTAINER_PROFILE="singularity"
else
  CONTAINER_PROFILE=""
fi

if [[ -n "$CONTAINER_PROFILE" ]]; then
  PROFILE_ARG="-profile ${NF_PROFILE},${CONTAINER_PROFILE}"
else
  echo "INFO: No Singularity/Apptainer detected — running without container profile"
  PROFILE_ARG="-profile ${NF_PROFILE}"
fi

############################################
# Run pipeline
############################################

cd "$(dirname "$0")/.." || die "Failed to change directory"

echo "Running: nextflow run $NF_PIPELINE $PROFILE_ARG"
nextflow run "$NF_PIPELINE" $PROFILE_ARG

############################################
# Locate latest test output
############################################

latest_testoutput=$(ls -1dt test_output* 2>/dev/null | head -n 1) \
  || die "No test_output directory found"

report_file=$(find "$latest_testoutput" -type f -name 'report*.csv' | head -n 1) \
  || die "No report CSV found"

############################################
# Compare outputs
############################################

generated_csv="$latest_testoutput/report_sorted.csv"
expected_csv="$latest_testoutput/report_ref_sorted.csv"

sort "$report_file" > "$generated_csv"
sort "$REF_CSV" > "$expected_csv"

if diff -q "$generated_csv" "$expected_csv" >/dev/null; then
  echo "PASSED ! CSV files match"
else
  echo "WARNING ! *** CSV files differ ***"
  diff "$generated_csv" "$expected_csv"
  exit 1
fi

############################################
# Cleanup
############################################

rm -f "$generated_csv" "$expected_csv"

# Fix terminal erase char after Nextflow
stty erase ^H || true

