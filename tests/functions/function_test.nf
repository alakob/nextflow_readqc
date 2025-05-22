#!/usr/bin/env nextflow

// Include the functions module
include { countReadsInFastq } from '../../functions.nf'

// Create test files directly in the workflow
workflow {
    // Create a test file with 2 reads in the current directory
    def testFile = file("test_reads.fastq")
    testFile.text = '''@read1
ACGTACGT
+
IIIIIIII
@read2
GCATGCAT
+
IIIIIIII
'''

    // Create an empty test file in the current directory
    def emptyFile = file("empty.fastq")
    emptyFile.text = ''

    // Test the function with both files
    def count1 = countReadsInFastq(testFile)
    def count2 = countReadsInFastq(emptyFile)
    
    // Print results
    println "Test file has ${count1} reads (expected: 2)"
    println "Empty file has ${count2} reads (expected: 0)"
    
    // Verify results
    assert count1 == 2, "Test file should have 2 reads but found ${count1}"
    assert count2 == 0, "Empty file should have 0 reads but found ${count2}"
    
    println "All tests passed successfully!"
}
