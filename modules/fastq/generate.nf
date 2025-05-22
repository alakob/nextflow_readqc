#!/usr/bin/env nextflow

// Process to generate synthetic FASTQ files using fastq-generator
process generate_fastq {
    // No container needed - using local fastq-generator installation

    output:
    tuple val("synthetic"), path("synthetic_R1.fastq"), path("synthetic_R2.fastq", optional: true)

    script:
    if (params.paired_end) {
        // For paired-end reads we need to run the command twice, once for each file
        """
        /Users/alakob/.local/bin/fastq-generator generate_random_fastq_se --sequence-size 100 --nb_seq ${params.num_reads ?: 1000} --output synthetic_R1.fastq
        /Users/alakob/.local/bin/fastq-generator generate_random_fastq_se --sequence-size 100 --nb_seq ${params.num_reads ?: 1000} --output synthetic_R2.fastq
        """
    } else {
        // For single-end, don't create R2 file
        """
        /Users/alakob/.local/bin/fastq-generator generate_random_fastq_se --sequence-size 100 --nb_seq ${params.num_reads ?: 1000} --output synthetic_R1.fastq
        touch synthetic_R2.fastq
        """
    }
}
