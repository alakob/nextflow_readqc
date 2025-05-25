#!/usr/bin/env nextflow

// Process to trim reads with fastp
process fastp_trim {
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(id), path(read1), path(read2)

    output:
    tuple val(id), path("trimmed_${id}_R1.fastq"), path("trimmed_${id}_R2.fastq", optional: true)

    script:
    if (read2) {
        """
        fastp -i ${read1} -I ${read2} -o trimmed_${id}_R1.fastq -O trimmed_${id}_R2.fastq --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --html /dev/null --json /dev/null
        """
    } else {
        """
        fastp -i ${read1} -o trimmed_${id}_R1.fastq --disable_quality_filtering --disable_length_filtering --disable_trim_poly_g --html /dev/null --json /dev/null
        """
    }
}
