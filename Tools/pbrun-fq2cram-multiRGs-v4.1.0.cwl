#!/usr/bin/env cwl-runner

class: CommandLineTool
id: pbrun-fq2cram-multiRGs
label: pbrun-fq2cram-multiRGs
cwlVersion: v1.1

$namespaces:
  cwltool: http://commonwl.org/cwltool#

requirements:
  DockerRequirement:
    dockerPull: hacchy/pbrun-fq2bam:4.1.0-1_v20231231
  ShellCommandRequirement: {}

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "11.0"
    cudaComputeCapability: "3.0"
    cudaDeviceCountMin: 1
    cudaDeviceCountMax: 8

baseCommand: [bash, /tools/pbrun-fq2bam-v4.1.0+.sh]

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

  knownSites:
    type: File[]?
    doc: A known indels file. The file must be in vcf.gz format. This option can be used multiple times.
    secondaryFiles:
      - .tbi
    inputBinding: 
      position: 10
      itemSeparator: ","

  low_memory:
    type: string?
    doc: Flag whether use --low-memory option
    default: "false"
    inputBinding:
      position: 9

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
  recal:
    type: File
    outputBinding:
      glob: $(inputs.prefix).bqsr.recal.table

arguments:
  - position: 8
    valueFrom: cram
