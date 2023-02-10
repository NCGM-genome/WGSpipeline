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
    type: File
    doc: FASTQ file 1

  fq2:
    type: File
    doc: FASTQ file 2

  RG_ID:
    type: string
    doc: Read group identifier (ID) in RG line
  RG_PL:
    type: string
    doc: Platform/technology used to produce the read (PL) in RG line
  RG_PU:
    type: string
    doc: Platform Unit (PU) in RG line
  RG_LB:
    type: string
    doc: DNA preparation library identifier (LB) in RG line
  RG_SM:
    type: string
    doc: Sample (SM) identifier in RG line  

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

  bwa_options:
    type: string?
    default: "-T 0 -Y"

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
    run: ../Tools/pbrun-fq2cram.cwl
    in:
      ref: ref
      fq1: fq1
      fq2: fq2
      RG_ID: RG_ID
      RG_PL: RG_PL
      RG_PU: RG_PU
      RG_LB: RG_LB
      RG_SM: RG_SM
      bwa_options: bwa_options
      num_gpus: num_gpus
      prefix: prefix
    out:
      - cram

  haplotypecaller_autosome:
    run: ../Tools/pbrun-haplotypecaller-from-cram.cwl
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
    run: ../Tools/pbrun-haplotypecaller-from-cram.cwl
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
    run: ../Tools/pbrun-haplotypecaller-from-cram.cwl
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
    run: ../Tools/pbrun-haplotypecaller-from-cram.cwl
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
    run: ../Tools/pbrun-haplotypecaller-from-cram.cwl
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
  
