#!/usr/bin/env nextflow
// reference https://seqera.io/training/#_get_started_with_nextflow


nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */

//params.baseDir = "/Users/sidd/Research/fun/nextflow"
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/transcriptome.fa"
params.outdir = "$baseDir/results"
params.multiqc = "$baseDir/multiqc"

println "reads: $params.reads"


log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         transcriptome: ${params.transcriptome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

process index {
    tag "indexing"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

process quantification {
    tag "quantifying" 
    publishDir "$baseDir", mode: 'copy'
    
    input:
    path index
    tuple val(pair_id), path(reads)
 
    output:
    path pair_id
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

process fastqc {
    tag "FASTQC on $sample_id"
    publishDir "$baseDir", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  

process multiqc {
    publishDir params.outdir, mode:'copy'
       
    input:
    path '*'
    
    output:
    path 'multiqc_report.html'
     
    script:
    """
    multiqc . 
    """
} 

workflow {
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    index(params.transcriptome)
    quantification(index.out, read_pairs_ch)
    fastqc (read_pairs_ch)
    multiqc(fastqc.out.collect())
}

workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}