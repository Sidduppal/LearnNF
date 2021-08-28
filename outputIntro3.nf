#!/usr/bin/env nextflow
// reference https://seqera.io/training/#_get_started_with_nextflow

nextflow.enable.dsl=2

process outputIntro3 {

    input:
    tuple val(sample_id), file(sample_files)

    output:
    tuple val(sample_id), file(sample_id)

    script:
    """
    echo your_command_here --reads $sample_id > $sample_id
    """
}

methods = channel.of('Hello world!', 'Yo, dude!', 'Duck!')

workflow {
    reads_ch = Channel.fromFilePairs('/Users/sidd/Research/fun/nextflow/data/ggal/*_{1,2}.fq')
    outputIntro3(reads_ch)
    outputIntro3.out.view{"Received: $it"}
}
