#!/usr/bin/env nextflow
// process run_GVCFtyper {
//   publishDir "${params.outdir}", mode:'copy'

//   input:
//   path GVCFLIST
//   path sh

//   output:
//   path ????

//   script:
//   """
//   $sh $GVCFLIST | perl -lnae "print $& if(/\d+/)"
//   """

// }
// process run_GVCFtyper_merge {
//   publishDir "${params.outdir}", mode:'copy'
// }
// process clean_space.sh {
//   publishDir "${params.outdir}", mode:'copy'
// }
// process run_varCal_SNP {
//   publishDir "${params.outdir}", mode:'copy'
// }
// process run_varCal_INDEL {
//   publishDir "${params.outdir}", mode:'copy'
// }
// process run_ApplyVarCal {
//   publishDir "${params.outdir}", mode:'copy'
// }

workflow {
  // params
  gvcfs_list = Channel.value(params.gvcfs_list)
  sh_run_GVCFtyper = Channel.value(params.sh_run_GVCFtyper)
  // shard to tuple
  shard = Channel
    .fromPath(params.shard)
    .splitText()
    .map { it.trim() }
    .toList()
    .flatMap { items -> 
        items.withIndex().collect { item, idx -> 
            def index = idx + 1
            def paddedIndex = index.toString().padLeft(3, '0')
            tuple(idx + 1, paddedIndex, item)
        }
    }
  shard.view()

  // run_GVCFtyper
  // out_run_GVCFtyper = run_GVCFtyper(gvcfs_list, sh_run_GVCFtyper)
}