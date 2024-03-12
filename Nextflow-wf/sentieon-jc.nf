#!/usr/bin/env nextflow
process run_GVCFtyper {
  publishDir "${params.outdir}", mode:'copy'

  input:
  path GVCFLIST
  path REF
  path REF_idx
  val SENTIEON_INSTALL_DIR
  val NCORE
  tuple val(idx), val(pad_idx), val(SHARD)

  output:
  tuple val(idx), val(pad_idx), path("*.vcf.gz")

  script:
  """
  ${SENTIEON_INSTALL_DIR}/bin/sentieon driver -t $NCORE -r $REF $SHARD \
    --traverse_param 2000/200 \
    --algo GVCFtyper \
    --genotype_mode multinomial \
    shard${pad_idx}.vcf.gz - < $GVCFLIST
  """
}

workflow {
  // params
  gvcfs_list = Channel.value(params.gvcfs_list)
  sh_GVCFtyper = Channel.value(params.run_GVCFtyper)
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

  ref = Channel.value(params.ref)
  sentieon = Channel.value(params.sentieon)
  ncore = Channel.value(params.ncore)
  ref_idx = Channel.value(params.ref_idx)

  // run_GVCFtyper
  out_GVCFtyper = run_GVCFtyper(gvcfs_list, ref, ref_idx, sentieon, ncore, shard)
}