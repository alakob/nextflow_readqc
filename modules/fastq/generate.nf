#!/usr/bin/env nextflow

// Process to generate synthetic FASTQ files using fastq-generator
process generate_fastq {
    // No container needed - using local fastq-generator installation

    output:
    tuple val("synthetic"), path("synthetic_R1.fastq"), path("synthetic_R2.fastq", optional: true)

    script:
    def num_reads = params.num_reads ?: 10  // Default to 10 reads instead of 1000
    def seq_size = params.seq_size ?: 100   // Make sequence size configurable
    def fastq_gen = params.fastq_generator ?: '/Users/alakob/.local/bin/fastq-generator'

    if (params.paired_end) {
        // For paired-end reads we need to run the command twice, once for each file
        """
        ${fastq_gen} generate_random_fastq_se --sequence-size ${seq_size} --nb_seq ${num_reads} --output synthetic_R1.fastq
        ${fastq_gen} generate_random_fastq_se --sequence-size ${seq_size} --nb_seq ${num_reads} --output synthetic_R2.fastq
        """
    } else {
        // For single-end, don't create R2 file
        """
        ${fastq_gen} generate_random_fastq_se --sequence-size ${seq_size} --nb_seq ${num_reads} --output synthetic_R1.fastq
        touch synthetic_R2.fastq
        """
    }
}
