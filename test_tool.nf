#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.reads = "/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz"

include{ FASTQC } from '/media/bigdrive1/sidd/autometa2_manuscript/metaBenchmarks/modules/nf-core/modules/fastqc/main.nf'

include {KRAKEN2_KRAKEN2} from '/media/bigdrive1/sidd/autometa2_manuscript/metaBenchmarks/modules/nf-core/modules/kraken2/kraken2/main.nf'

workflow {
    reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    db=file("/media/bigdrive1/Databases/kraken2/kraken2_db")
    FASTQC(reads)
    KRAKEN2_KRAKEN2 (reads, db)
}