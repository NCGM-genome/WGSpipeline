#!/usr/bin/env nextflow
params.ncore = 2

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
process cnvkit_batch {
  publishDir params.outdir, mode:'copy'

  input:
  path out_bam
  path out_bai
  path ref_cnn
  val sample
  val ncore
  
  output:
  path "${sample}-diagram.pdf", emit: pdf
  path "${sample}-scatter.png", emit: png
  path "${sample}.antitargetcoverage.cnn", emit: anti
  path "${sample}.targetcoverage.cnn", emit: target
  path "${sample}.bintest.cns", emit: cns_bin
  path "${sample}.call.cns", emit: cns_call
  path "${sample}.cns", emit: cns
  path "${sample}.cnr", emit: cnr

  script:
  """
  cnvkit.py \\
        batch \\
        $out_bam \\
        -r $ref_cnn \\
        --scatter --diagram \\
        -m wgs --segment-method cbs \\
        -p $ncore \\
  """
   // -d $outdir --scatter --diagram \\
}

workflow {
  // params
  ref = Channel.fromPath(params.ref)
  cram = Channel.fromPath(params.cram)
  sample = Channel.of(params.sample)
  ref_cnn = Channel.fromPath(params.ref_cnn)
  ncore = Channel.of(params.ncore)
  // cram2bam
  cram2bam_out =cram2bam(ref, cram, sample)
  out_bam = cram2bam_out.bam
  // bam2bai
  bam2bai_out = bam2bai(out_bam)
  out_bai = bam2bai_out.bai
  // batch
  batch_out = cnvkit_batch(out_bam, out_bai, ref_cnn, sample ,ncore)
}
