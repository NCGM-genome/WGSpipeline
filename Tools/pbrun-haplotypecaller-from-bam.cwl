#!/usr/bin/env cwl-runner

class: CommandLineTool
id: pbrun-haplotypecaller
label: pbrun-haplotypecaller
cwlVersion: v1.1

$namespaces:
  cwltool: http://commonwl.org/cwltool#

requirements:
  DockerRequirement:
    dockerPull: nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1
  ShellCommandRequirement: {}

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "11.0"
    cudaComputeCapabilityMin: "3.0"
    deviceCountMin: 1
    deviceCountMax: 8

baseCommand: [pbrun, haplotypecaller]

inputs:
  bam:
    type: File
    doc: BAM file
    inputBinding:
      prefix: --in-bam
      position: 2
    secondaryFiles:
      - .bai

  ref: 
    type: File
    doc: Reference FASTA file
    inputBinding:
      prefix: --ref
      position: 1
    secondaryFiles:
      - ^.dict
      - .fai

  interval_file:
    type: File
    inputBinding:
      prefix: --interval-file
      position: 3

  ploidy:
    type: int
    inputBinding:
      prefix: --ploidy
      position: 4

  num_gpus:
    type: int
    inputBinding:
      prefix: --num-gpus
      position: 5

  prefix:
    type: string
    doc: Output file prefix

outputs:
  gvcf:
    type: File
    outputBinding:
      glob: $(inputs.prefix).g.vcf.gz
    secondaryFiles:
      - .tbi

arguments:
  - position: 6
    prefix: --out-variants
    valueFrom: $(inputs.prefix).g.vcf.gz
  - position: 7
    valueFrom: --gvcf

