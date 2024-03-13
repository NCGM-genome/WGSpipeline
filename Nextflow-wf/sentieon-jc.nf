#!/usr/bin/env nextflow
process run_GVCFtyper {
  publishDir "${params.outdir}/shard", mode:'copy'

  input:
  path GVCFLIST
  path REF
  path REF_idx
  val NCORE
  tuple val(idx), val(pad_idx), val(SHARD)

  output:
  tuple val(idx), val(pad_idx), path("*.vcf.gz"), path("*.vcf.gz.tbi")

  script:
  """
  ${params.sentieon}/bin/sentieon driver -t $NCORE -r $REF $SHARD \
    --traverse_param 2000/200 \
    --algo GVCFtyper \
    --genotype_mode multinomial \
    shard${pad_idx}.vcf.gz - < $GVCFLIST
  """
}

// process run_GVCFtyper_merge {
//   publishDir "${params.outdir}", mode:'copy'

//   input:
//   path REF
//   path REF_idx
//   val NCORE
//   tuple val(cn), val(vcfs)

//   output:

//   script:
//   """
//   ${params.sentieon}/bin/sentieon driver --passthru -t $NCORE -r $REF \
//     --algo GVCFtyper --merge \
//     tmp.GVCFtyper.chr${cn}.vcf.gz $vcfs
//   """
// }

workflow {
  // GVCFtyper
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
  ncore = Channel.value(params.ncore)
  ref_idx = Channel.value(params.ref_idx)
  // run
  out_GVCFtyper = run_GVCFtyper(gvcfs_list, ref, ref_idx, ncore, shard)


  // GVCFtyper_merge
  // params
  // shard_chr = Channel
  //   .fromPath(params.shard_chr)
  //   .splitText()
  //   .map { it.trim() }
  //   .flatMap { line ->
  //       // 空白で分割して各数値を処理
  //       def paddedNumbers = line.split(/\s+/).collect { num ->
  //           // 数値を3桁のゼロパディングでフォーマットし、前後に文字列を追加
  //           "${params.outdir}/shard/shard${num.toString().padLeft(3, '0')}.vcf.gz"
  //       }.join(' ') // パディングを適用した数値を空白で結合
  //   }    
  //   .toList()
  //   .flatMap { items -> 
  //       items.withIndex().collect { item, idx -> 
  //           def index = idx + 1
  //           tuple(index, item)
  //       }
  //   }
  //   .flatMap { tuple ->
  //       def index = tuple[0] // 元の1要素目の値を取得
  //       tuple[1].split(' ').collect { file ->
  //           // 新しいタプルを生成
  //           [index, file]
  //       }
  //   }
 
  // // run
  // shard_chr.view()

   table = Channel
    .fromPath(params.shard_chr)
    .splitText()
    .map { it.trim() }   
    .toList()
    .flatMap { items -> 
        items.withIndex().collect { item, idx -> 
            def index = idx + 1
            tuple(index.toString(), item)
        }
    }
    .flatMap { tuple ->
        def index = tuple[0] // 元の1要素目の値を取得
        tuple[1].split(' ').collect { file ->
            // 新しいタプルを生成
            [file, index]
        }
    } 
  table.view()
  // out_GVCFtyper.view()
  paired = table.join(out_GVCFtyper, by: [0, 0])
  paired.view()
  // out_GVCFtyper_merge = run_GVCFtyper_merge(ref, ref_idx, ncore, shard_chr)
}