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

process configManta {
  publishDir "${params.outdir}/${out_bam.baseName}", mode:'copy'

  input:
  path out_bam
  path out_bai
  path region
  path region_idx
  path ref
  path ref_idx
  val ncore
  
  output:
  // path "${params.outdir}/${out_bam.baseName}/runWorkflow.py", emit: runWorkflow
  path "MantaWorkflow/runWorkflow.py", emit: runWorkflow
  path "*.vcf", emit: vcf

  script:
  """
  configManta.py \
    --bam $out_bam \
    --referenceFasta $ref \
    --callRegions $region
  MantaWorkflow/runWorkflow.py -j $ncore
  """
  //     --runDir ${params.outdir}/${out_bam.baseName} \
}

// process runWorkflow {
//   publishDir "${params.outdir}/${sample}", mode:'copy'

//   input:
//   // path out_py
//   val sample
//   val ncore
  
//   output:
//   path "*.vcf"

//   script:
//   """
//   ${params.outdir}/${sample}/MantaWorkflow/runWorkflow.py -j $ncore
//   """
//   // ${params.outdir}/${sample}/MantaWorkflow/runWorkflow.py -j $ncore
// }

workflow {
  // params -> Channel
  sample = Channel
    .fromPath(params.sample_list)
    .splitText()
    .map { it.trim() }
    .map { file ->
      def fileName = file.split('/').last().replace(".cram", "")
      return fileName
    }
  cram = Channel
    .fromPath(params.sample_list)
    .splitText()
    .map { it.trim() }
  ref = Channel.value(params.ref)
  ref_idx = Channel.value(params.ref_idx)
  region = Channel.value(params.region)
  region_idx = Channel.value(params.region_idx)
  ncore = Channel.value(params.ncore)

  // cram2bam process
  cram2bam_out =cram2bam(ref, cram, sample)
  out_bam = cram2bam_out.bam
  // bam2bai process
  bam2bai_out = bam2bai(out_bam)
  out_bai = bam2bai_out.bai

  // configManta process
  configManta_out = configManta(out_bam, out_bai, region, region_idx ,ref ,ref_idx, ncore)
  // out_py = configManta_out.runWorkflow

  // // runWorkflow process
  // runWorkflow_out = runWorkflow(out_py, ncore)
  // runWorkflow_out = runWorkflow(sample, ncore)
}
