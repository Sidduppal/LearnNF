#!/usr/bin/env nextflow
// reference https://seqera.io/training/#_get_started_with_nextflow

nextflow.enable.dsl=2

process outputIntro {

    input:
    val x

    output:
    val x

    script:
    """
    echo $x
    """
}

methods = channel.of('Hello world!', 'Yo, dude!', 'Duck!')

workflow {
    outputIntro(methods)
    outputIntro.out.view{"Received: $it"}
}
