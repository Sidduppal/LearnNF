#!/usr/bin/env nextflow
// reference https://seqera.io/training/#_get_started_with_nextflow
nextflow.enable.dsl=2

process concatenateFile {
    echo true 

    input:
    file reads

    output:
    stdout

    script:

    """
    cat $reads > concat_file.fq
    head -n 20 concat_file.fq
    """
}

reads = Channel.fromPath( '/Users/sidd/Research/fun/nextflow/data/ggal/*.fq' )

workflow {

    concatenateFile(reads.collect())

}