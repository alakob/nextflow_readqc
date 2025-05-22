# Nextflow Tutorial Pipeline

A demonstration Nextflow pipeline for processing FASTQ files, performing quality control, and trimming.

## Features

- Generate synthetic FASTQ files or use existing ones
- Process both paired-end and single-end reads
- Run FastQC on raw reads
- Trim reads with fastp
- Run FastQC on trimmed reads
- Count reads in files using a custom function

## Pipeline Workflows

The pipeline includes several workflows:

- **Main workflow**: Autodetects inputs and processes them accordingly
- **generate_single_end**: Generates and processes single-end reads
- **generate_paired_end**: Generates and processes paired-end reads
- **existing_single_end**: Processes existing single-end reads
- **existing_paired_end**: Processes existing paired-end reads

## Requirements

- Nextflow 22.10.0 or later
- Docker (optional, but recommended)
- FastQC
- fastp
- fastq-generator (for synthetic read generation)

## Usage

```bash
# Run with default parameters (generate paired-end reads)
nextflow run main.nf

# Use existing paired-end reads
nextflow run main.nf --generate_reads false --paired_end true --reads "path/to/reads/*_{1,2}.fastq"

# Use existing single-end reads
nextflow run main.nf --generate_reads false --paired_end false --single_end_reads "path/to/reads/*.fastq"
```

## Testing

This pipeline includes comprehensive tests using the nf-test framework:

```bash
# Run all tests
nf-test test

# Run specific test
nf-test test tests/processes/fastqc.nf.test
```

## CI/CD

The pipeline includes GitHub Actions workflows for continuous integration and testing.

## License

MIT
