#!/usr/bin/env nextflow

// Process to run FastQC on a single read file
process fastqc {
    tag "$id:$read"
    // Container directive moved to nextflow.config
    publishDir path: { params.publish_dir ?: publish_dir }, mode: 'copy'
    
    input:
    val publish_dir
    tuple val(id), path(read)
    
    output:
    path "${read.baseName}_fastqc.zip"
    
    script:
    """
    fastqc ${read}
    """
}
