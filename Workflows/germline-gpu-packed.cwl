{
    "$graph": [
        {
            "class": "CommandLineTool",
            "id": "#pbrun-fq2cram-multiRGs.cwl",
            "label": "pbrun-fq2cram-multiRGs",
            "requirements": [
                {
                    "dockerPull": "hacchy/pbrun-fq2bam:4.0.0-1_v20230412",
                    "class": "DockerRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "cudaVersionMin": "11.0",
                    "cudaComputeCapability": "3.0",
                    "cudaDeviceCountMin": 1,
                    "cudaDeviceCountMax": 8,
                    "class": "http://commonwl.org/cwltool#CUDARequirement"
                }
            ],
            "baseCommand": [
                "bash",
                "/tools/pbrun-fq2bam.sh"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "-T 0 -Y",
                    "inputBinding": {
                        "position": 5
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/bwa_options"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "FASTQ file 1",
                    "inputBinding": {
                        "position": 1,
                        "itemSeparator": ","
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/fq1"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "FASTQ file 2",
                    "inputBinding": {
                        "position": 2,
                        "itemSeparator": ","
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/fq2"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "A known indels file. The file must be in vcf.gz format. This option can be used multiple times.",
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "inputBinding": {
                        "position": 9,
                        "itemSeparator": ","
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/knownSites"
                },
                {
                    "type": "int",
                    "inputBinding": {
                        "position": 6
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/num_gpus"
                },
                {
                    "type": "string",
                    "doc": "Output file prefix",
                    "inputBinding": {
                        "position": 7
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/prefix"
                },
                {
                    "type": "File",
                    "doc": "Reference FASTA file",
                    "secondaryFiles": [
                        {
                            "pattern": "^.dict",
                            "required": null
                        },
                        {
                            "pattern": ".fai",
                            "required": null
                        },
                        {
                            "pattern": ".64.amb",
                            "required": null
                        },
                        {
                            "pattern": ".64.ann",
                            "required": null
                        },
                        {
                            "pattern": ".64.bwt",
                            "required": null
                        },
                        {
                            "pattern": ".64.pac",
                            "required": null
                        },
                        {
                            "pattern": ".64.sa",
                            "required": null
                        },
                        {
                            "pattern": ".64.alt",
                            "required": null
                        }
                    ],
                    "inputBinding": {
                        "position": 4
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/ref"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "doc": "Read group string",
                    "inputBinding": {
                        "position": 3,
                        "itemSeparator": ","
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/rg"
                }
            ],
            "arguments": [
                {
                    "position": 8,
                    "valueFrom": "cram"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.prefix).cram"
                    },
                    "secondaryFiles": [
                        {
                            "pattern": ".crai",
                            "required": null
                        }
                    ],
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/cram"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.prefix).bqsr.recal.table"
                    },
                    "id": "#pbrun-fq2cram-multiRGs.cwl/pbrun-fq2cram-multiRGs/recal"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "id": "#pbrun-haplotypecaller-from-cram.cwl",
            "label": "pbrun-haplotypecaller",
            "requirements": [
                {
                    "dockerPull": "nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1",
                    "class": "DockerRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "hints": [
                {
                    "cudaVersionMin": "11.0",
                    "cudaComputeCapability": "3.0",
                    "cudaDeviceCountMin": 1,
                    "cudaDeviceCountMax": 8,
                    "class": "http://commonwl.org/cwltool#CUDARequirement"
                }
            ],
            "baseCommand": [
                "pbrun",
                "haplotypecaller"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "BAM file",
                    "inputBinding": {
                        "prefix": "--in-bam",
                        "position": 2
                    },
                    "secondaryFiles": [
                        {
                            "pattern": ".crai",
                            "required": null
                        }
                    ],
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/cram"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--interval-file",
                        "position": 3
                    },
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/interval_file"
                },
                {
                    "type": "int",
                    "inputBinding": {
                        "prefix": "--num-gpus",
                        "position": 5
                    },
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/num_gpus"
                },
                {
                    "type": "int",
                    "inputBinding": {
                        "prefix": "--ploidy",
                        "position": 4
                    },
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/ploidy"
                },
                {
                    "type": "string",
                    "doc": "Output file prefix",
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/prefix"
                },
                {
                    "type": "File",
                    "doc": "Reference FASTA file",
                    "inputBinding": {
                        "prefix": "--ref",
                        "position": 1
                    },
                    "secondaryFiles": [
                        {
                            "pattern": "^.dict",
                            "required": null
                        },
                        {
                            "pattern": ".fai",
                            "required": null
                        }
                    ],
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/ref"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.prefix).g.vcf.gz"
                    },
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "id": "#pbrun-haplotypecaller-from-cram.cwl/pbrun-haplotypecaller/gvcf"
                }
            ],
            "arguments": [
                {
                    "position": 6,
                    "prefix": "--out-variants",
                    "valueFrom": "$(inputs.prefix).g.vcf.gz"
                },
                {
                    "position": 7,
                    "valueFrom": "--gvcf"
                }
            ]
        },
        {
            "class": "Workflow",
            "id": "#main",
            "label": "germline",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Interval BED file for PAR regions",
                    "id": "#main/PAR_interval"
                },
                {
                    "type": "File",
                    "doc": "Interval BED file for autosome regions",
                    "id": "#main/autosome_interval"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "default": "-T 0 -Y",
                    "id": "#main/bwa_options"
                },
                {
                    "type": "File",
                    "doc": "Interval BED file for chrX regions",
                    "id": "#main/chrX_interval"
                },
                {
                    "type": "File",
                    "doc": "Interval BED file for chrY regions",
                    "id": "#main/chrY_interval"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "FASTQ file 1",
                    "id": "#main/fq1"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "FASTQ file 2",
                    "id": "#main/fq2"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "A known indels file. The file must be in vcf.gz format. This option can be used multiple times.",
                    "secondaryFiles": [
                        {
                            "pattern": ".tbi",
                            "required": null
                        }
                    ],
                    "id": "#main/knownSites"
                },
                {
                    "type": "int",
                    "id": "#main/num_gpus"
                },
                {
                    "type": "string",
                    "doc": "Output file prefix",
                    "id": "#main/prefix"
                },
                {
                    "type": "File",
                    "doc": "Reference FASTA file",
                    "secondaryFiles": [
                        {
                            "pattern": "^.dict",
                            "required": null
                        },
                        {
                            "pattern": ".fai",
                            "required": null
                        },
                        {
                            "pattern": ".64.amb",
                            "required": null
                        },
                        {
                            "pattern": ".64.ann",
                            "required": null
                        },
                        {
                            "pattern": ".64.bwt",
                            "required": null
                        },
                        {
                            "pattern": ".64.pac",
                            "required": null
                        },
                        {
                            "pattern": ".64.sa",
                            "required": null
                        },
                        {
                            "pattern": ".64.alt",
                            "required": null
                        }
                    ],
                    "id": "#main/ref"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "doc": "Read group string",
                    "id": "#main/rg"
                }
            ],
            "steps": [
                {
                    "run": "#pbrun-fq2cram-multiRGs.cwl",
                    "in": [
                        {
                            "source": "#main/bwa_options",
                            "id": "#main/fq2cram/bwa_options"
                        },
                        {
                            "source": "#main/fq1",
                            "id": "#main/fq2cram/fq1"
                        },
                        {
                            "source": "#main/fq2",
                            "id": "#main/fq2cram/fq2"
                        },
                        {
                            "source": "#main/knownSites",
                            "id": "#main/fq2cram/knownSites"
                        },
                        {
                            "source": "#main/num_gpus",
                            "id": "#main/fq2cram/num_gpus"
                        },
                        {
                            "source": "#main/prefix",
                            "id": "#main/fq2cram/prefix"
                        },
                        {
                            "source": "#main/ref",
                            "id": "#main/fq2cram/ref"
                        },
                        {
                            "source": "#main/rg",
                            "id": "#main/fq2cram/rg"
                        }
                    ],
                    "out": [
                        "#main/fq2cram/cram",
                        "#main/fq2cram/recal"
                    ],
                    "id": "#main/fq2cram"
                },
                {
                    "run": "#pbrun-haplotypecaller-from-cram.cwl",
                    "in": [
                        {
                            "source": "#main/fq2cram/cram",
                            "id": "#main/haplotypecaller_PAR/cram"
                        },
                        {
                            "source": "#main/PAR_interval",
                            "id": "#main/haplotypecaller_PAR/interval_file"
                        },
                        {
                            "source": "#main/num_gpus",
                            "id": "#main/haplotypecaller_PAR/num_gpus"
                        },
                        {
                            "valueFrom": "$(2)",
                            "id": "#main/haplotypecaller_PAR/ploidy"
                        },
                        {
                            "valueFrom": "$(inputs.tmpprefix).PAR",
                            "id": "#main/haplotypecaller_PAR/prefix"
                        },
                        {
                            "source": "#main/ref",
                            "id": "#main/haplotypecaller_PAR/ref"
                        },
                        {
                            "source": "#main/prefix",
                            "id": "#main/haplotypecaller_PAR/tmpprefix"
                        }
                    ],
                    "out": [
                        "#main/haplotypecaller_PAR/gvcf"
                    ],
                    "id": "#main/haplotypecaller_PAR"
                },
                {
                    "run": "#pbrun-haplotypecaller-from-cram.cwl",
                    "in": [
                        {
                            "source": "#main/fq2cram/cram",
                            "id": "#main/haplotypecaller_autosome/cram"
                        },
                        {
                            "source": "#main/autosome_interval",
                            "id": "#main/haplotypecaller_autosome/interval_file"
                        },
                        {
                            "source": "#main/num_gpus",
                            "id": "#main/haplotypecaller_autosome/num_gpus"
                        },
                        {
                            "valueFrom": "$(2)",
                            "id": "#main/haplotypecaller_autosome/ploidy"
                        },
                        {
                            "valueFrom": "$(inputs.tmpprefix).autosome",
                            "id": "#main/haplotypecaller_autosome/prefix"
                        },
                        {
                            "source": "#main/ref",
                            "id": "#main/haplotypecaller_autosome/ref"
                        },
                        {
                            "source": "#main/prefix",
                            "id": "#main/haplotypecaller_autosome/tmpprefix"
                        }
                    ],
                    "out": [
                        "#main/haplotypecaller_autosome/gvcf"
                    ],
                    "id": "#main/haplotypecaller_autosome"
                },
                {
                    "run": "#pbrun-haplotypecaller-from-cram.cwl",
                    "in": [
                        {
                            "source": "#main/fq2cram/cram",
                            "id": "#main/haplotypecaller_chrX_female/cram"
                        },
                        {
                            "source": "#main/chrX_interval",
                            "id": "#main/haplotypecaller_chrX_female/interval_file"
                        },
                        {
                            "source": "#main/num_gpus",
                            "id": "#main/haplotypecaller_chrX_female/num_gpus"
                        },
                        {
                            "valueFrom": "$(2)",
                            "id": "#main/haplotypecaller_chrX_female/ploidy"
                        },
                        {
                            "valueFrom": "$(inputs.tmpprefix).chrX_female",
                            "id": "#main/haplotypecaller_chrX_female/prefix"
                        },
                        {
                            "source": "#main/ref",
                            "id": "#main/haplotypecaller_chrX_female/ref"
                        },
                        {
                            "source": "#main/prefix",
                            "id": "#main/haplotypecaller_chrX_female/tmpprefix"
                        }
                    ],
                    "out": [
                        "#main/haplotypecaller_chrX_female/gvcf"
                    ],
                    "id": "#main/haplotypecaller_chrX_female"
                },
                {
                    "run": "#pbrun-haplotypecaller-from-cram.cwl",
                    "in": [
                        {
                            "source": "#main/fq2cram/cram",
                            "id": "#main/haplotypecaller_chrX_male/cram"
                        },
                        {
                            "source": "#main/chrX_interval",
                            "id": "#main/haplotypecaller_chrX_male/interval_file"
                        },
                        {
                            "source": "#main/num_gpus",
                            "id": "#main/haplotypecaller_chrX_male/num_gpus"
                        },
                        {
                            "valueFrom": "$(1)",
                            "id": "#main/haplotypecaller_chrX_male/ploidy"
                        },
                        {
                            "valueFrom": "$(inputs.tmpprefix).chrX_male",
                            "id": "#main/haplotypecaller_chrX_male/prefix"
                        },
                        {
                            "source": "#main/ref",
                            "id": "#main/haplotypecaller_chrX_male/ref"
                        },
                        {
                            "source": "#main/prefix",
                            "id": "#main/haplotypecaller_chrX_male/tmpprefix"
                        }
                    ],
                    "out": [
                        "#main/haplotypecaller_chrX_male/gvcf"
                    ],
                    "id": "#main/haplotypecaller_chrX_male"
                },
                {
                    "run": "#pbrun-haplotypecaller-from-cram.cwl",
                    "in": [
                        {
                            "source": "#main/fq2cram/cram",
                            "id": "#main/haplotypecaller_chrY/cram"
                        },
                        {
                            "source": "#main/chrY_interval",
                            "id": "#main/haplotypecaller_chrY/interval_file"
                        },
                        {
                            "source": "#main/num_gpus",
                            "id": "#main/haplotypecaller_chrY/num_gpus"
                        },
                        {
                            "valueFrom": "$(1)",
                            "id": "#main/haplotypecaller_chrY/ploidy"
                        },
                        {
                            "valueFrom": "$(inputs.tmpprefix).chrY",
                            "id": "#main/haplotypecaller_chrY/prefix"
                        },
                        {
                            "source": "#main/ref",
                            "id": "#main/haplotypecaller_chrY/ref"
                        },
                        {
                            "source": "#main/prefix",
                            "id": "#main/haplotypecaller_chrY/tmpprefix"
                        }
                    ],
                    "out": [
                        "#main/haplotypecaller_chrY/gvcf"
                    ],
                    "id": "#main/haplotypecaller_chrY"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputSource": "#main/fq2cram/cram",
                    "id": "#main/cram"
                },
                {
                    "type": "File",
                    "outputSource": "#main/haplotypecaller_PAR/gvcf",
                    "id": "#main/gvcf_PAR"
                },
                {
                    "type": "File",
                    "outputSource": "#main/haplotypecaller_autosome/gvcf",
                    "id": "#main/gvcf_autosome"
                },
                {
                    "type": "File",
                    "outputSource": "#main/haplotypecaller_chrX_female/gvcf",
                    "id": "#main/gvcf_chrX_female"
                },
                {
                    "type": "File",
                    "outputSource": "#main/haplotypecaller_chrX_male/gvcf",
                    "id": "#main/gvcf_chrX_male"
                },
                {
                    "type": "File",
                    "outputSource": "#main/haplotypecaller_chrY/gvcf",
                    "id": "#main/gvcf_chrY"
                },
                {
                    "type": "File",
                    "outputSource": "#main/fq2cram/recal",
                    "id": "#main/recal"
                }
            ]
        }
    ],
    "cwlVersion": "v1.1"
}
