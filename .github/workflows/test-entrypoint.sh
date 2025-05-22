#!/bin/bash
set -e

# Function to print section headers
section() {
    echo "===================================================================="
    echo "= $1"
    echo "===================================================================="
}

# Navigate to mounted pipeline directory
cd /pipeline

# Print diagnostic information
section "Environment Info"
echo "Nextflow version: $(nextflow -v)"
echo "Python version: $(python --version)"
echo "nf-test version: $(nf-test --version)"
echo "fastq-generator version: $(fastq-generator --version || echo 'Version command not supported')"
echo "FastQC version: $(fastqc --version)"

# Run function tests directly with Nextflow
section "Running Function Tests"
nextflow run tests/functions/function_test.nf

# Run nf-test suite
section "Running nf-test Suite"
nf-test test --profile=test

echo "All tests completed!"
