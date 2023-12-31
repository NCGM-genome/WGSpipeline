#!/usr/bin/env cwl-runner

class: Workflow
id: germline
label: germline
cwlVersion: v1.1

requirements:
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  fq1:
    type: File[]
    doc: FASTQ file 1

  fq2:
    type: File[]
    doc: FASTQ file 2

  rg:
    type: string[]
    doc: Read group string

  ref: 
    type: File
    doc: Reference FASTA file
    secondaryFiles:
      - ^.dict
      - .fai
      - .64.amb
      - .64.ann
      - .64.bwt
      - .64.pac
      - .64.sa
      - .64.alt

  knownSites:
    type: File[]?
    doc: A known indels file. The file must be in vcf.gz format. This option can be used multiple times.
    secondaryFiles:
      - .tbi

  low_memory:
    type: string?
    doc: Flag whether use --low-memory option
    default: "false"

  bwa_options:
    type: string?
    default: "-Y"

  num_gpus:
    type: int

  prefix:
    type: string
    doc: Output file prefix

  autosome_interval:
    type: File
    doc: Interval BED file for autosome regions

  PAR_interval:
    type: File
    doc: Interval BED file for PAR regions

  chrX_interval:
    type: File
    doc: Interval BED file for chrX regions

  chrY_interval:
    type: File
    doc: Interval BED file for chrY regions

steps:
  fq2cram:
    run: ../Tools/pbrun-fq2cramfast-multiRGs-v4.2.1.cwl
    in:
      ref: ref
      knownSites: knownSites
      fq1: fq1
      fq2: fq2
      rg: rg
      low_memory: low_memory
      bwa_options: bwa_options
      num_gpus: num_gpus
      prefix: prefix
    out:
      - cram
      - recal

  haplotypecaller_autosome:
    run: ../Tools/pbrun-haplotypecaller-from-cram-v4.2.1.cwl
    in: 
      ref: ref
      cram: fq2cram/cram
      interval_file: autosome_interval
      ploidy: 
        valueFrom: $(2)
      num_gpus: num_gpus
      tmpprefix: prefix
      prefix: 
        valueFrom: $(inputs.tmpprefix).autosome
    out: 
      - gvcf

  haplotypecaller_PAR:
    run: ../Tools/pbrun-haplotypecaller-from-cram-v4.2.1.cwl
    in: 
      ref: ref
      cram: fq2cram/cram
      interval_file: PAR_interval
      ploidy: 
        valueFrom: $(2)
      num_gpus: num_gpus
      tmpprefix: prefix
      prefix: 
        valueFrom: $(inputs.tmpprefix).PAR
    out: 
      - gvcf

  haplotypecaller_chrX_female:
    run: ../Tools/pbrun-haplotypecaller-from-cram-v4.2.1.cwl
    in: 
      ref: ref
      cram: fq2cram/cram
      interval_file: chrX_interval
      ploidy: 
        valueFrom: $(2)
      num_gpus: num_gpus
      tmpprefix: prefix
      prefix: 
        valueFrom: $(inputs.tmpprefix).chrX_female
    out: 
      - gvcf

  haplotypecaller_chrX_male:
    run: ../Tools/pbrun-haplotypecaller-from-cram-v4.2.1.cwl
    in: 
      ref: ref
      cram: fq2cram/cram
      interval_file: chrX_interval
      ploidy: 
        valueFrom: $(1)
      num_gpus: num_gpus
      tmpprefix: prefix
      prefix: 
        valueFrom: $(inputs.tmpprefix).chrX_male
    out: 
      - gvcf

  haplotypecaller_chrY:
    run: ../Tools/pbrun-haplotypecaller-from-cram-v4.2.1.cwl
    in: 
      ref: ref
      cram: fq2cram/cram
      interval_file: chrY_interval
      ploidy: 
        valueFrom: $(1)
      num_gpus: num_gpus
      tmpprefix: prefix
      prefix: 
        valueFrom: $(inputs.tmpprefix).chrY
    out: 
      - gvcf

outputs:
  cram:
    type: File
    outputSource: fq2cram/cram

  recal:
    type: File
    outputSource: fq2cram/recal

  gvcf_autosome:
    type: File
    outputSource: haplotypecaller_autosome/gvcf

  gvcf_PAR:
    type: File
    outputSource: haplotypecaller_PAR/gvcf

  gvcf_chrX_female:
    type: File
    outputSource: haplotypecaller_chrX_female/gvcf

  gvcf_chrX_male:
    type: File
    outputSource: haplotypecaller_chrX_male/gvcf

  gvcf_chrY:
    type: File
    outputSource: haplotypecaller_chrY/gvcf
  
