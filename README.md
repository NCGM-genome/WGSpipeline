# WGSpipeline

This repository contains the following workflow:
- **`germline-gpu.cwl`:** This workflow calculates `cram` and `gvcf` files from `fastq` files using `parabricks` version 4.0.0.

## Installation requirements
- [Hardware requirements to run parabricks](https://docs.nvidia.com/clara/parabricks/4.0.0/GettingStarted.html#hardware-requirements)
- [Software requirements to run parabricks](https://docs.nvidia.com/clara/parabricks/4.0.0/GettingStarted.html#software-requirements)
- SingularityCE 3.10.2+
- python 3+

## Install `cwltool`
You can install cwltool by executing the following commands:
```
$ cd /path/to/working/directory/
$ python -m venv cwlenv
$ . cwlenv/bin/activate
$ python -m pip install --upgrade pip
$ pip install cwltool
```

Please confirm the cwltool version:
```
$ cwltool --version
/path/to/working/directory/cwlenv/bin/cwltool 3.1.20230201224320
```


## Usage

### `germline-gpu.cwl` workflow
```
usage: Workflows/germline-gpu.cwl [-h] [--bwa_options BWA_OPTIONS] \
                                    --ref REF \
                                    --fq1 FQ1 \
                                    --fq2 FQ2 \
                                    --rg RG \
                                    --autosome_interval AUTOSOME_INTERVAL \
                                    --PAR_interval PAR_INTERVAL \
                                    --chrX_interval CHRX_INTERVAL \
                                    --chrY_interval CHRY_INTERVAL \
                                    --num_gpus NUM_GPUS \
                                    --prefix PREFIX 

optional arguments:
  -h, --help                Show this help message and exit.
  --bwa_options STRING      Pass supported bwa mem options as one string. 
                            The current original bwa mem supported options are -M, -Y, and -T. 
                            (e.g. --bwa-options="-M -Y") 
                            (default: "-T 0 -Y")
  --ref FILE                Path to the reference file.
  --fq1 FILE                Path to FASTQ file 1.
  --fq2 FILE                Path to FASTQ file 2.
  --rg STRING               Read group string.
  --autosome_interval FILE  Path to interval BED file for autosome regions.
  --PAR_interval FILE       Path to interval BED file for PAR regions.
  --chrX_interval FILE      Path to interval BED file for chrX regions.
  --chrY_interval FILE      Path to interval BED file for chrY regions.
  --num_gpus INT            Number of GPUs to use for a run (should be ≥1). 
  --prefix STRING           Output file prefix.
```

Basic usage of `germline-gpu.cwl` workflow with a pair of FASTQ files.
```
$ cd /path/to/working/directory
$ . cwlenv/bin/activate
$ cwltool --singularity \
    --outdir output_directory \
    WGSpipeline/Workflows/germline-gpu.cwl \
    --ref reference_hg38/Homo_sapiens_assembly38.fasta \
    --fq1 libraryA_1.fq.gz \
    --fq2 libraryA_2.fq.gz \
    --rg "@RG\\tID:read_group_id\\tPL:platform\\tPU:platform_unit\\tLB:library\\tSM:sample_id" \
    --num_gpus 4 \
    --prefix prefix \
    --autosome_interval WGSpipeline/interval_files/autosome.bed \
    --PAR_interval WGSpipeline/interval_files/PAR.bed \
    --chrX_interval WGSpipeline/interval_files/chrX.bed \
    --chrY_interval WGSpipeline/interval_files/chrY.bed
```

Usage of `germline-gpu.cwl` workflow with multiple pairs of FASTQ files.
```
$ cd /path/to/working/directory
$ . cwlenv/bin/activate
$ cwltool --singularity \
    --outdir output_directory \
    WGSpipeline/Workflows/germline-gpu.cwl \
    --ref reference_hg38/Homo_sapiens_assembly38.fasta \
    --fq1 libraryA_1.fq.gz \
    --fq2 libraryA_2.fq.gz \
    --rg "@RG\\tID:read_group_id\\tPL:platform\\tPU:platform_unit\\tLB:libraryA\\tSM:sample_id" \
    --fq1 libraryB_1.fq.gz \
    --fq2 libraryB_2.fq.gz \
    --rg "@RG\\tID:read_group_id\\tPL:platform\\tPU:platform_unit\\tLB:libraryB\\tSM:sample_id" \
    --num_gpus 4 \
    --prefix prefix \
    --autosome_interval WGSpipeline/interval_files/autosome.bed \
    --PAR_interval WGSpipeline/interval_files/PAR.bed \
    --chrX_interval WGSpipeline/interval_files/chrX.bed \
    --chrY_interval WGSpipeline/interval_files/chrY.bed
```
`--fq1`, `--fq2`, and `--rg` options can be repeated multiple times. 
The number of `--fq1` options should be same as the number of `--fq2` and `--rg` options.
The order of `--fq1`, `--fq2`, and `--rg` options are important. 



