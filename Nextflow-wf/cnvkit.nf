#!/usr/bin/env nextflow
process cram2bam {
  publishDir "${params.outdir}/${sample}", mode:'copy'
  label 'samtools'

  input:
  path ref
  val cram
  val sample
  
  output:
  path "${sample}.bam", emit: bam

  script:
  """
  samtools view -@ ${params.samtools_thread} -b -T $ref $cram > ${sample}.bam
  """
}

process bam2bai {
  publishDir "${params.outdir}/${out_bam.baseName}", mode:'copy'
  label 'samtools'

  input:
  path out_bam
  
  output:
  path "${out_bam}.bai", emit: bai

  script:
  """
  samtools index $out_bam
  """
}

process cnvkit_batch {
  publishDir "${params.outdir}/${out_bam.baseName}", mode:'copy'

  input:
  path out_bam
  path out_bai
  path cnn
  val ncore
  
  output:
  path "*.pdf", emit: pdf
  path "*.png", emit: png
  path "*.antitargetcoverage.cnn", emit: anti
  path "*.targetcoverage.cnn", emit: target
  path "*.bintest.cns", emit: cns_bin
  path "*.call.cns", emit: cns_call
  path "*.cns", emit: cns
  path "*.cnr", emit: cnr

  script:
  """
  cnvkit.py \\
        batch \\
        $out_bam \\
        -r $cnn \\
        --scatter --diagram \\
        -m wgs --segment-method cbs \\
        -p $ncore \\
  """
}

workflow {
  // params
  sample = Channel.fromPath(params.sample_list)
                      .splitText()
                      .map { it.trim() }
                      .map { it.replace(".cram", "") }
  cram = Channel.fromPath(params.sample_list)
                      .splitText()
                      .map { it.trim() }
                      .map { "${params.base_path}${it}" }  
  ref = Channel.value(params.ref)
  cnn = Channel.value(params.cnn)
  ncore = Channel.value(params.ncore)

  // cram2bam process
  cram2bam_out =cram2bam(ref, cram, sample)
  out_bam = cram2bam_out.bam
  // bam2bai process
  bam2bai_out = bam2bai(out_bam)
  out_bai = bam2bai_out.bai
  // batch process
  batch_out = cnvkit_batch(out_bam, out_bai, cnn ,ncore)
}
