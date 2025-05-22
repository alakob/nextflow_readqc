#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.reads = "data/reads/*_{1,2}.fastq"  // Path to existing reads (if not generating)
params.paired_end = true                      // Toggle for paired-end or single-end reads
params.outdir = "results"                     // Output directory
params.generate_reads = true                // Set to false to use existing reads instead of generating
params.single_end_reads = null                // Added for existing_single_end workflow

// Import modules
include { generate_fastq } from './modules/fastq/generate.nf'
include { fastp_trim } from './modules/fastp/trim.nf'
include { fastqc_wf as fastqc_raw_wf; fastqc_wf as fastqc_trimmed_wf } from './modules/fastqc/workflow.nf'

// Import utility functions
include { countReadsInFastq } from './functions.nf'

// Main workflow definition
workflow {
    // Choose input: generate synthetic reads or use provided ones
    if (params.generate_reads) {
        generate_fastq()
        input_channel = generate_fastq.out
    } else {
        if (params.paired_end) {
            input_channel = Channel.fromFilePairs(params.reads, size: 2).map { id, files -> tuple(id, files[0], files[1]) }
        } else {
            input_channel = Channel.fromPath(params.reads).map { file -> tuple(file.baseName, file, null) }
        }
    }

    // Process the data using the shared sub-workflow
    process_reads(input_channel)
}

// Define a shared sub-workflow for the common processing steps
workflow process_reads {
    take:
        input_channel  // Input channel of reads (id, read1, read2)
    
    main:
        // Count reads in input files using our utility function
        input_channel.subscribe { id, r1, r2 ->
            def r1_count = countReadsInFastq(r1)
            println "[READ COUNT] Input file ${r1.name}: ${r1_count} reads"
            
            if (r2) {
                def r2_count = countReadsInFastq(r2)
                println "[READ COUNT] Input file ${r2.name}: ${r2_count} reads"
            }
        }
        
        // Process the data - this was duplicated across all workflows
        input_reads = input_channel.flatMap { id, r1, r2 -> r2 ? [ [id, r1], [id, r2] ] : [ [id, r1] ] }
        fastqc_raw_wf("${params.outdir}/fastqc_raw", input_reads)
        
        trimming = fastp_trim(input_channel)
        
        // Count reads in trimmed files
        trimming.subscribe { id, tr1, tr2 ->
            def tr1_count = countReadsInFastq(tr1)
            println "[READ COUNT] Trimmed file ${tr1.name}: ${tr1_count} reads"
            
            if (tr2) {
                def tr2_count = countReadsInFastq(tr2)
                println "[READ COUNT] Trimmed file ${tr2.name}: ${tr2_count} reads"
            }
        }
        
        trimmed_reads = trimming.flatMap { id, tr1, tr2 -> tr2 ? [ [id, tr1], [id, tr2] ] : [ [id, tr1] ] }
        fastqc_trimmed_wf("${params.outdir}/fastqc_trimmed", trimmed_reads)
    
    emit:
        trimmed = trimming
}

// Define named workflows for testing
workflow generate_single_end {
    // Set parameters
    params.paired_end = false
    params.generate_reads = true
    
    // Generate reads
    generate_fastq()
    input_channel = generate_fastq.out
    
    // Process the data using the shared sub-workflow
    process_reads(input_channel)
}

workflow generate_paired_end {
    // Set parameters
    params.paired_end = true
    params.generate_reads = true
    
    // Generate reads
    generate_fastq()
    input_channel = generate_fastq.out
    
    // Process the data using the shared sub-workflow
    process_reads(input_channel)
}

workflow existing_single_end {
    // Set parameters
    params.paired_end = false
    
    // Override generate_reads for testing purposes if needed
    // This allows the test to use synthetic reads while still testing this workflow
    if (params.generate_reads) {
        // Use the generate_fastq process to create synthetic reads
        generate_fastq()
        input_channel = generate_fastq.out
    } else {
        // Get existing reads - use single_end_reads for single-end mode
        input_channel = Channel.fromPath(params.single_end_reads ?: params.reads).map { singleEndFile -> 
            // Use an empty file for the second read
            def dummyFile = file('empty.fastq')
            if (!dummyFile.exists()) {
                dummyFile.text = ''
            }
            [singleEndFile.simpleName, singleEndFile, dummyFile]
        }
    }
    
    // Process the data using the shared sub-workflow
    process_reads(input_channel)
}

workflow existing_paired_end {
    // Set parameters
    params.paired_end = true
    params.generate_reads = false
    
    // Get existing reads
    input_channel = Channel.fromFilePairs(params.reads, size: 2).map { id, files -> tuple(id, files[0], files[1]) }
    
    // Process the data using the shared sub-workflow
    process_reads(input_channel)
}

