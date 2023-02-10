#!/usr/bin/env cwl-runner

class: CommandLineTool
id: pbrun-fq2cram-multiRGs
label: pbrun-fq2cram-multiRGs
cwlVersion: v1.1

$namespaces:
  cwltool: http://commonwl.org/cwltool#

requirements:
  DockerRequirement:
    dockerPull: hacchy/pbrun-fq2bam:4.0.0-1_v20230210
  ShellCommandRequirement: {}

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "11.4"
    cudaComputeCapabilityMin: "3.0"
    deviceCountMin: 1
    deviceCountMax: 8

baseCommand: [bash, /tools/pbrun-fq2bam.sh]

inputs:
  fq1:
    type: File[]
    doc: FASTQ file 1
    inputBinding:
      position: 1
      itemSeparator: ","

  fq2:
    type: File[]
    doc: FASTQ file 2
    inputBinding:
      position: 2
      itemSeparator: ","

  rg:
    type: string[]
    doc: Read group string
    inputBinding:
      position: 3
      itemSeparator: ","

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
    inputBinding:
      position: 4

  bwa_options:
    type: string?
    default: "-T 0 -Y"
    inputBinding: 
      position: 5

  num_gpus:
    type: int
    inputBinding:
      position: 6

  prefix:
    type: string
    doc: Output file prefix
    inputBinding:
      position: 7

outputs:
  cram:
    type: File
    outputBinding:
      glob: $(inputs.prefix).cram
    secondaryFiles:
      - .crai

arguments:
  - position: 8
    valueFrom: cram
