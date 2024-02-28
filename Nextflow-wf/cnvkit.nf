#!/usr/bin/env nextflow
process cram2bam {
  publishDir params.outdir, mode:'copy'
  label 'samtools'

  input:
    path ref
    path cram
    val sample
  
  output:
  path "${sample}.bam", emit: bam

  script:
  """
  samtools view -b -T $ref $cram > ${sample}.bam
  """
}
process bam2bai {
  publishDir params.outdir, mode:'copy'
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

params.ncore = '2'
process batch {
  publishDir params.outdir, mode:'copy'
  label 'cnvkit'

  input:
    path out_bam
  
  output:
  path "${out_bam}.bai", emit: bai

  script:
  """
  cnvkit.py \\
        batch \\
        $bam \\
        $fasta_args \\
        $reference_args \\
        --processes $task.cpus \\
        $args
  """
}

workflow {
  // params
  ref = Channel.fromPath(params.ref)
  cram = Channel.fromPath(params.cram)
  sample = Channel.of(params.sample)
  // cram2bam
  cram2bam_out =cram2bam(ref, cram, sample)
  out_bam = cram2bam_out.bam
  // bam2bai
  bam2bai_out = bam2bai(out_bam)
  out_bai = bam2bai_out.bai
}
