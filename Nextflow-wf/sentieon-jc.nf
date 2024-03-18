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
  tuple val(cn), path(vcf), path("${params.prefix}.chr${cn}.vcf.gz.tbi")

  script:
  """
  ${params.sentieon}/bin/sentieon util vcfindex ${vcf} && \
  rm ${params.outdir}/tmp.GVCFtyper.chr${cn}.vcf.gz ${params.outdir}/tmp.GVCFtyper.chr${cn}.vcf.gz.tbi
  """
}

process merge_siteonly {
  publishDir "${params.outdir}", mode:'copy'
  label 'bamtols'

  input:
  tuple val(cn), path(vcf), path(vcf_tbi)

  output:
  tuple val(cn), path("${params.prefix}.siteonly*vcf.gz")

  script:
  """
  if [ $cn -eq 1 ]; then
      bcftools view --threads ${params.ncore} -Oz -G ${params.prefix}.chr${cn}.vcf.gz > ${params.prefix}.siteonly.vcf.gz
  else
      bcftools view --threads ${params.ncore} -Oz -GH ${params.prefix}.chr${cn}.vcf.gz > ${params.prefix}.siteonly.chr${cn}.vcf.gz
  fi
  """
}

process clean_space {
  publishDir "${params.outdir}", mode:'copy'
  label 'tabix'

  input:
  path siteonly_vcfs

  output:
  path "${params.prefix}.siteonly.vcf.gz", emit: vcf
  path "${params.prefix}.siteonly.vcf.gz.tbi", emit: idx

  script:
  """
  for cn in {2..22}
  do
      cat ${params.prefix}.siteonly.chr\${cn}.vcf.gz >> ${params.prefix}.siteonly.vcf.gz
      rm ${params.prefix}.siteonly.chr\${cn}.vcf.gz
      rm ${params.outdir}/${params.prefix}.siteonly.chr\${cn}.vcf.gz
  done
  rm -Rf ${params.outdir}/shard
  tabix ${params.prefix}.siteonly.vcf.gz
  """
}

process varCal_SNP {
  publishDir "${params.outdir}", mode:'copy'

  input:
  path REF
  path ref_idx
  path SITEONLY
  path SITEONLY_idx
  path DBSNP
  path DBSNP_idx
  path HAPMAP
  path HAPMAP_idx
  path OMNI
  path OMNI_idx
  path i1000G
  path i1000G_idx
  
  output:
  path "*.tranches", emit: tranches
  path "*.R", emit: R
  path "*.recal", emit: recal
  path "*.recal.idx", emit: recal_idx

  script:
  """
  mode="SNP"
  ANN="--annotation QD --annotation MQRankSum --annotation ReadPosRankSum --annotation FS --annotation SOR --annotation DP"
  MAXGAUSSIAN=6
  RESOURCE="
    --resource $HAPMAP --resource_param hapmap,known=false,training=true,truth=true,prior=15 \
    --resource $OMNI --resource_param omni,known=false,training=true,truth=true,prior=12 \
    --resource $i1000G --resource_param 1000G,known=false,training=true,truth=false,prior=10 \
    --resource $DBSNP --resource_param dbsnp,known=true,training=false,truth=false,prior=7"
  TRANCHE="--tranche 100.0 --tranche 99.95 --tranche 99.9 --tranche 99.8 --tranche 99.6 --tranche 99.5 --tranche 99.4 --tranche 99.3 --tranche 99.0 --tranche 98.0 --tranche 97.0 --tranche 90.0"
  ${params.sentieon}/bin/sentieon driver -t ${params.ncore} -r $REF \
    --algo VarCal \
    -v $SITEONLY \
    --max_gaussians \${MAXGAUSSIAN} \
    --var_type \${mode} \
    \${TRANCHE} \${ANN} \${RESOURCE} \
    --tranches_file \${mode}.tranches \
    --plot_file VarCal.\${mode}.R \
    VarCal.\${mode}.recal
  """
}

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
  ref_idx = Channel.value(params.ref_idx)
  //// run
  out_GVCFtyper = run_GVCFtyper(gvcfs_list, ref, ref_idx, shard)

  // GVCFtyper_merge
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
  //// run
  out_GVCFtyper_merge = run_GVCFtyper_merge(ref, ref_idx, shard_chr, out_GVCFtyper.collect())
  //// merge_bcftools
  out_merge_bcftools = merge_bcftools(out_GVCFtyper_merge)
  //// merge_sentieon_idx
  out_merge_sentieon_idx = merge_sentieon_idx(out_merge_bcftools)
  //// merge_siteonly
  out_merge_siteonly = merge_siteonly(out_merge_sentieon_idx)

  // clean_space
  siteonly_vcfs = out_merge_siteonly.map { it[1] }.collect()
  out_clean_space = clean_space(siteonly_vcfs)
  siteonly_vcf = out_clean_space.vcf
  siteonly_idx = out_clean_space.idx

  // varCal-SNP
  dbsnp_path =Channel.fromPath("${params.bundle}/${params.dbsnp}")
  dbsnp_path_idx =Channel.fromPath("${params.bundle}/${params.dbsnp}.tbi")

  hapmap_path =Channel.fromPath("${params.bundle}/${params.hapmap}")
  hapmap_path_idx =Channel.fromPath("${params.bundle}/${params.hapmap}.tbi")

  omni_path =Channel.fromPath("${params.bundle}/${params.omni}")
  omni_path_idx =Channel.fromPath("${params.bundle}/${params.omni}.tbi")

  i1000g_path =Channel.fromPath("${params.bundle}/${params.i1000g}")
  i1000g_path_idx =Channel.fromPath("${params.bundle}/${params.i1000g}.tbi")

  //// run
  out_varCal_SNP = varCal_SNP(ref, 
                              ref_idx, 
                              siteonly_vcf, 
                              siteonly_idx, 
                              dbsnp_path, 
                              dbsnp_path_idx, 
                              hapmap_path, 
                              hapmap_path_idx, 
                              omni_path, 
                              omni_path_idx, 
                              i1000g_path,
                              i1000g_path_idx)
}
