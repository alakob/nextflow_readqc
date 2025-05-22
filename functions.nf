#!/usr/bin/env nextflow

// Simple function to count sequence lines in a FASTQ file
// Returns the number of sequences (reads) in the file
def countReadsInFastq(fastqFile) {
    // For testing, handle strings by converting them to Files
    def file = fastqFile instanceof String ? file(fastqFile) : fastqFile
    
    // Guard against non-existent files
    if (!file.exists()) {
        return 0
    }
    
    def lines = file.readLines()
    return lines.size() / 4 // Each FASTQ record is 4 lines
}

// Return functions
return [
    countReadsInFastq: this.&countReadsInFastq
]
