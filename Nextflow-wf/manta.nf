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

// process configManta {
//   publishDir "${params.outdir}/${sample_name}", mode:'copy'

//   input:
//   path region
//   path region_idx
//   path ref
//   path ref_idx
//   tuple val(sample_name), path(bam_path), path(bai_path)
  
//   output:
//   tuple val(sample_name), path("MantaWorkflow/*.py")

//   script:
//   """
//   configManta.py \
//     --bam $bam_path \
//     --referenceFasta $ref \
//     --callRegions $region
  
//   """
// }

process configManta {
  publishDir "${params.outdir}/${sample_name}", mode:'copy'

  input:
  path region
  path region_idx
  path ref
  path ref_idx
  val ncore
  tuple val(sample_name), path(bam_path), path(bai_path)
  
  output:
  // tuple val(sample_name), path("MantaWorkflow/*.py")
  path "MantaWorkflow/**"

  script:
  """
  configManta.py \
    --bam $bam_path \
    --referenceFasta $ref \
    --callRegions $region
  MantaWorkflow/runWorkflow.py -j $ncore
  """
}

// process runWorkflow {
//   publishDir "${params.outdir}/${sample_name}", mode:'copy'

//   input:
//   val ncore
//   tuple val(sample_name), path(py_path)
  
//   output:
//   // path "MantaWorkflow/results/variants*.vcf"
//   path "MantaWorkflow/**"

//   script:
//   """
//   MantaWorkflow/$py_path -j $ncore
//   """
// }

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
  region = Channel.value(params.region)
  region_idx = Channel.value(params.region_idx)
  ncore = Channel.value(params.ncore)


  // cram2bam process
  cram2bam_out =cram2bam(ref, ref_idx, cram)
  cram2bam_out.view()
  
  // bam2bai process
  bam2bai_out = bam2bai(cram2bam_out)
  bam2bai_out.view()
  // tuple join
  paired = cram2bam_out.join(bam2bai_out, by: 0)
  paired.view { "${it[0]}, ${it[1]}, ${it[2]}" }

  // configManta process
  // configManta_out = configManta(out_bam, out_bai, region, region_idx ,ref ,ref_idx, ncore)
  // configManta_out = configManta(region, region_idx ,ref ,ref_idx, paired)
  configManta_out = configManta(region, region_idx ,ref ,ref_idx, ncore, paired)
  configManta_out.view()

  // // runWorkflow process
  // runWorkflow_out = runWorkflow(ncore, configManta_out)
  // runWorkflow_out.view()
  
}
