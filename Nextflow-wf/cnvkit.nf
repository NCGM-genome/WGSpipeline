#!/usr/bin/env nextflow
process cram2bam {
  publishDir "${params.outdir}/${sample_name}", mode:'copy'
  label 'samtools'

  input:
  path ref
  path ref_idx
  tuple val(sample_name), path(cram_path)

  output:
  tuple val(sample_name), path("*.bam")

  script:
  """
  samtools view -@ ${params.samtools_thread} -b -T $ref $cram_path > ${sample_name}.bam
  """
}

process bam2bai {
  publishDir "${params.outdir}/${sample_name}", mode:'copy'
  label 'samtools'

  input:
  tuple val(sample_name), path(bam_path)
  
  output:
  tuple val(sample_name), path("*.bai")

  script:
  """
  samtools index $bam_path
  """
}

process cnvkit_batch {
  publishDir "${params.outdir}/${sample_name}", mode:'copy'

  input:
  path cnn
  val ncore
  tuple val(sample_name), path(bam_path), path(bai_path)
  
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
        $bam_path \\
        -r $cnn \\
        --scatter --diagram \\
        -m wgs --segment-method cbs \\
        -p $ncore \\
  """
}

workflow {
  // params -> Channel
  cram = Channel
    .fromPath(params.sample_list)
    .splitText()
    .map { it.trim() }
    .map {
        def sampleName = it.split('/').last().replace(".cram", "")
        tuple(sampleName, it)
    } 
  ref = Channel.value(params.ref)
  ref_idx = Channel.value(params.ref_idx)
  cnn = Channel.value(params.cnn)
  ncore = Channel.value(params.ncore)

  // cram2bam process
  cram2bam_out =cram2bam(ref, ref_idx, cram)
  // bam2bai process
  bam2bai_out = bam2bai(cram2bam_out)
  // tuple join
  paired = cram2bam_out.join(bam2bai_out, by: 0)
  // batch process
  batch_out = cnvkit_batch(cnn ,ncore, paired)
}

// workflow {
//   // params -> Channel
//   sample = Channel
//     .fromPath(params.sample_list)
//     .splitText()
//     .map { it.trim() }
//     .map { file ->
//       def fileName = file.split('/').last().replace(".cram", "")
//       return fileName
//     }
//   cram = Channel
//     .fromPath(params.sample_list)
//     .splitText()
//     .map { it.trim() }
//   ref = Channel.value(params.ref)
//   cnn = Channel.value(params.cnn)
//   ncore = Channel.value(params.ncore)

//   // cram2bam process
//   cram2bam_out =cram2bam(ref, cram, sample)
//   out_bam = cram2bam_out.bam
//   // bam2bai process
//   bam2bai_out = bam2bai(out_bam)
//   out_bai = bam2bai_out.bai
//   // batch process
//   batch_out = cnvkit_batch(out_bam, out_bai, cnn ,ncore)
// }
