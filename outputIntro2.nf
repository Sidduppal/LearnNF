#!/usr/bin/env nextflow
// reference https://seqera.io/training/#_get_started_with_nextflow

nextflow.enable.dsl=2

process outputIntro2 {

    output:
    file 'result.txt'

    script:
    """
    echo \$RANDOM > result.txt
    """
}

workflow {
    outputIntro2()
    outputIntro2.out.view{"Received: $it"}
}