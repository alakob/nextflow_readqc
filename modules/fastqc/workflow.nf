#!/usr/bin/env nextflow

include { fastqc } from './run.nf'

// Define a workflow for running FastQC with a specific output directory
workflow fastqc_wf {
    take:
        publish_dir    // Directory for output
        reads_channel  // Channel of reads
        
    main:
        // Run the process with the specified output directory
        fastqc(publish_dir, reads_channel)
        
    emit:
        fastqc.out
}
