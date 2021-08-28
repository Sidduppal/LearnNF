#!/usr/bin/env nextflow
// reference https://seqera.io/training/#_get_started_with_nextflow

nextflow.enable.dsl=2

process foo {
  echo true
  
  input:
  val x  
  val y 
  
  script:
   """
   echo $x and $y
   """
}

num = Channel.from(1,2,3)
let = Channel.from('a','b','c')

workflow {

    foo(num, let)

}