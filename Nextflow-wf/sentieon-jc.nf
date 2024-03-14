#!/usr/bin/env nextflow
process run_GVCFtyper {
  publishDir "${params.outdir}/shard", mode:'copy'

  input:
  path GVCFLIST
  path REF
  path REF_idx
  tuple val(idx), val(pad_idx), val(SHARD)

  output:
  tuple val(idx), path("*.vcf.gz"), path("*.vcf.gz.tbi")

  script:
  """
  ${params.sentieon}/bin/sentieon driver -t ${params.ncore}  -r $REF $SHARD \
    --traverse_param 2000/200 \
    --algo GVCFtyper \
    --genotype_mode multinomial \
    shard${pad_idx}.vcf.gz - < $GVCFLIST
  """
}

process run_GVCFtyper_merge {
  publishDir "${params.outdir}", mode:'copy'

  input:
  path REF
  path REF_idx
  tuple val(cn), val(vcfs)
  val tuple_collects

  output:
  tuple val(cn), path("tmp.GVCFtyper.chr${cn}.vcf.gz"), path("tmp.GVCFtyper.chr${cn}.vcf.gz.tbi") 

  script:
  """
  ${params.sentieon}/bin/sentieon driver --passthru -t ${params.ncore}  -r $REF \
    --algo GVCFtyper --merge \
    tmp.GVCFtyper.chr${cn}.vcf.gz $vcfs
  """
}

process merge_bcftools {
  publishDir "${params.outdir}", mode:'copy'
  label 'bamtols'

  input:
  tuple val(cn), path(tmp_vcf), path(tmp_vcf_tbi)

  output:
  tuple val(cn), path("${params.prefix}.chr${cn}.vcf.gz")

  script:
  """
  bcftools view --threads ${params.ncore} -r chr${cn} -Oz -o ${params.prefix}.chr${cn}.vcf.gz tmp.GVCFtyper.chr${cn}.vcf.gz
  """
}

process merge_sentieon_idx {
  publishDir "${params.outdir}", mode:'copy'

  input:
  tuple val(cn), path(vcf)

  output:
  tuple val(cn), path("${params.prefix}.chr${cn}.vcf.gz.tbi")

  script:
  """
  ${params.sentieon}/bin/sentieon util vcfindex ${vcf} && \
  rm ${params.outdir}/tmp.GVCFtyper.chr${cn}.vcf.gz ${params.outdir}/tmp.GVCFtyper.chr${cn}.vcf.gz.tbi
  """
}

workflow {
  //// GVCFtyper
  gvcfs_list = Channel.value(params.gvcfs_list)
  sh_GVCFtyper = Channel.value(params.run_GVCFtyper)
  shard = Channel
    .fromPath(params.shard)
    .splitText()
    .map { it.trim() }
    .toList()
    .flatMap { items -> 
        items.withIndex().collect { item, idx -> 
            def index = idx + 1
            def paddedIndex = index.toString().padLeft(3, '0')
            tuple(index, paddedIndex, item)
        }
    }
  ref = Channel.value(params.ref)
  ref_idx = Channel.value(params.ref_idx)
  // run
  out_GVCFtyper = run_GVCFtyper(gvcfs_list, ref, ref_idx, shard)

  //// GVCFtyper_merge
  shard_chr = Channel
    .fromPath(params.shard_chr)
    .splitText()
    .map { it.trim() }
    .flatMap { line ->
        // 空白で分割して各数値を処理
        def paddedNumbers = line.split(/\s+/).collect { num ->
            // 数値を3桁のゼロパディングでフォーマットし、前後に文字列を追加
            "${params.outdir}/shard/shard${num.toString().padLeft(3, '0')}.vcf.gz"
        }.join(' ') // パディングを適用した数値を空白で結合
    }    
    .toList()
    .flatMap { items -> 
        items.withIndex().collect { item, idx -> 
            def index = idx + 1
            tuple(index, item)
        }
    }
  // run
  out_GVCFtyper_merge = run_GVCFtyper_merge(ref, ref_idx, shard_chr, out_GVCFtyper.collect())
  // merge_bcftools
  out_merge_bcftools = merge_bcftools(out_GVCFtyper_merge)

  // merge_sentieon_idx
  out_merge_sentieon_idx = merge_sentieon_idx(out_merge_bcftools)
}
