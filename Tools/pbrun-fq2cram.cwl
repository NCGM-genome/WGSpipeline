#!/usr/bin/env cwl-runner

class: CommandLineTool
id: pbrun-fq2cram
label: pbrun-fq2cram
cwlVersion: v1.1

$namespaces:
  cwltool: http://commonwl.org/cwltool#

requirements:
  DockerRequirement:
    dockerPull: nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1
  ShellCommandRequirement: {}

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "11.4"
    cudaComputeCapabilityMin: "3.0"
    deviceCountMin: 1
    deviceCountMax: 8

baseCommand: [pbrun, fq2bam]

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
    inputBinding:
      prefix: --ref
      position: 1
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
    inputBinding:
      prefix: --bwa-options
      position: 4
      shellQuote: true

  num_gpus:
    type: int
    inputBinding:
      prefix: --num-gpus
      position: 5

  prefix:
    type: string
    doc: Output file prefix

outputs:
  cram:
    type: File
    outputBinding:
      glob: $(inputs.prefix).cram
    secondaryFiles:
      - .crai

arguments:
  - position: 6
    prefix: --out-bam
    valueFrom: $(inputs.prefix).cram
  - position: 7
    prefix: --in-fq
    valueFrom: $(inputs.fq1.path) $(inputs.fq2.path) "@RG\\tID:$(inputs.RG_ID)\\tLB:$(inputs.RG_LB)\\tPL:$(inputs.RG_PL)\\tSM:$(inputs.RG_SM)\\tPU:$(inputs.RG_PU)"
    shellQuote: false

