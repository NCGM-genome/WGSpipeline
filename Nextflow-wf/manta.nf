#!/usr/bin/env nextflow
process cram2bam {
  publishDir "${params.outdir}/${cram.baseName}", mode:'copy'
  label 'samtools'

  input:
  path ref
  path cram
  
  output:
  // path "*.bam", emit: bam
  tuple val(cram.baseName), path("*.bam"), emit: bam

  script:
  def baseName = file(cram).getBaseName()
  """
  samtools view -@ ${params.samtools_thread} -b -T $ref $cram > ${baseName}.bam
  """
}

process bam2bai {
  publishDir "${params.outdir}/${out_bam.baseName}", mode:'copy'
  label 'samtools'

  input:
  path out_bam
  
  output:
  path "*.bai", emit: bai

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

// workflow {
//   // params -> Channel
//   cram = Channel
//     .fromPath(params.sample_list)
//     .splitText()
//     .map { it.trim() }
//   ref = Channel.value(params.ref)
//   ref_idx = Channel.value(params.ref_idx)
//   region = Channel.value(params.region)
//   region_idx = Channel.value(params.region_idx)
//   ncore = Channel.value(params.ncore)

//   // cram2bam process
//   cram2bam_out =cram2bam(ref, cram)
//   out_bam = cram2bam_out.bam
//   out_bam.view()
//   // bam2bai process
//   bam2bai_out = bam2bai(out_bam)
//   out_bai = bam2bai_out.bai
//    out_bai.view()

//   // // configManta process
//   // configManta_out = configManta(out_bam, out_bai, region, region_idx ,ref ,ref_idx, ncore)
//   // out_py = configManta_out.runWorkflow

//   // // runWorkflow process
//   // runWorkflow_out = runWorkflow(out_py, ncore)
//   // runWorkflow_out = runWorkflow(sample, ncore)
// }

workflow {
  // Parameters -> Channel
  Channel
    .fromPath(params.sample_list)
    .splitText()
    .map { it.trim() }
    .set { cram_files }

  ref = file(params.ref)

  // Process cram2bam and bam2bai in a single flow
  cram_files
    .map { cram -> tuple(ref, cram) }
    .flatMap { ref, cram ->
      cram2bam(ref, cram)
        .flatMap { bam ->
          bam2bai(bam)
            .map { bai -> tuple(bam, bai) }
        }
    }
    .view { tuple ->
      "BAM file: ${tuple[0]}, BAI file: ${tuple[1]}"
    }
}
