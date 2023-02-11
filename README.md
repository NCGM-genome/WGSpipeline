# WGSpipeline

## Usage

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



