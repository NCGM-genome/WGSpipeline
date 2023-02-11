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

## Clone this repository
Clone this repository by executing the following commands:
```
$ cd /path/to/working/directory/
$ git clone https://github.com/NCGM-genome/WGSpipeline.git
```

## Download reference and resource files
Download reference and resource files from the URLs listed in [reference_hg38.download_links.txt](./download_links/reference_hg38.download_links.txt) by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=reference_hg38 ; mkdir -p $OUTDIR ; for url in `cat WGSpipeline/download_links/reference_hg38.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```

## Usage of `germline-gpu.cwl` workflow
```
usage: Workflows/germline-gpu.cwl [-h] [--bwa_options STRING] \
                                    --ref FILE \
				    [--knownSites FILE] \
                                    --fq1 FILE \
                                    --fq2 FILE \
                                    --rg STRING \
                                    --autosome_interval FILE \
                                    --PAR_interval FILE \
                                    --chrX_interval FILE \
                                    --chrY_interval FILE \
                                    --num_gpus INT \
                                    --prefix STRING

optional arguments:
  -h, --help                Show this help message and exit.
  --bwa_options STRING      Pass supported bwa mem options as one string. 
                            The current original bwa mem supported options are -M, -Y, and -T. 
                            (e.g. --bwa-options="-M -Y") 
                            (default: "-T 0 -Y")
  --ref FILE                Path to the reference file.	
  --knownSites FILE         Path to a known indels file. 
  	       		    The file must be in vcf.gz format. 
			    This option can be used multiple times.
  --fq1 FILE                Path to FASTQ file 1.
			    This option can be used multiple times.
  --fq2 FILE                Path to FASTQ file 2.
			    This option can be used multiple times.
  --rg STRING               Read group string.
			    This option can be used multiple times.
  --autosome_interval FILE  Path to interval BED file for autosome regions.
  --PAR_interval FILE       Path to interval BED file for PAR regions.
  --chrX_interval FILE      Path to interval BED file for chrX regions.
  --chrY_interval FILE      Path to interval BED file for chrY regions.
  --num_gpus INT            Number of GPUs to use for a run (should be ≥1). 
  --prefix STRING           Output file prefix.
```

## Download a dataset for tutorial
Download a dataset for tutorial from the URLs listed in [wgs_fastq_NA12878_20k.download_links.txt](./download_links/wgs_fastq_NA12878_20k.download_links.txt) by executing the following commnds:
```
$ cd /path/to/working/directory/
$ OUTDIR=wgs_fastq ; mkdir -p $OUTDIR ; for url in `cat WGSpipeline/download_links/wgs_fastq_NA12878_20k.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```


## Tutorial 1: Run  `germline-gpu.cwl` workflow with a pair of FASTQ files.
Run `germline-gpu.cwl` workflow with a pair of FASTQ files.

```
$ cd /path/to/working/directory
$ mkdir -p tutorial_01
$ . cwlenv/bin/activate
$ cwltool --singularity \
    --outdir tutorial_01 \
    WGSpipeline/Workflows/germline-gpu.cwl \
    --ref reference_hg38/Homo_sapiens_assembly38.fasta \
    --fq1 wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06HDADXX130110.1\\tPL:ILLUMINA\\tPU:H06HDADXX130110.1\\tLB:H06HDADXX130110.1\\tSM:NA12878" \
    --num_gpus 4 \
    --prefix NA12878.H06HDADXX130110.1 \
    --autosome_interval WGSpipeline/interval_files/autosome.bed \
    --PAR_interval WGSpipeline/interval_files/PAR.bed \
    --chrX_interval WGSpipeline/interval_files/chrX.bed \
    --chrY_interval WGSpipeline/interval_files/chrY.bed
```
When `--knownSites` option is not used, then BQSR will not be applied and `_prefix_.bqsr.recla.table` file will be empty. 

Output files will be saved in the directory `/path/to/working/directory/tutorial_01`:
```
/path/to/working/directory/tutorial_01
|--NA12878.H06HDADXX130110.1.PAR.g.vcf.gz
|--NA12878.H06HDADXX130110.1.PAR.g.vcf.gz.tbi
|--NA12878.H06HDADXX130110.1.autosome.g.vcf.gz
|--NA12878.H06HDADXX130110.1.autosome.g.vcf.gz.tbi
|--NA12878.H06HDADXX130110.1.bqsr.recal.table
|--NA12878.H06HDADXX130110.1.chrX_female.g.vcf.gz
|--NA12878.H06HDADXX130110.1.chrX_female.g.vcf.gz.tbi
|--NA12878.H06HDADXX130110.1.chrX_male.g.vcf.gz
|--NA12878.H06HDADXX130110.1.chrX_male.g.vcf.gz.tbi
|--NA12878.H06HDADXX130110.1.chrY.g.vcf.gz
|--NA12878.H06HDADXX130110.1.chrY.g.vcf.gz.tbi
|--NA12878.H06HDADXX130110.1.cram
|--NA12878.H06HDADXX130110.1.cram.crai
```


## Tutorial 2: Run `germline-gpu.cwl` workflow with multiple pairs of FASTQ files.
Run `germline-gpu.cwl` workflow with multiple pairs of FASTQ files.

`--fq1`, `--fq2`, and `--rg` options can be repeated multiple times. 
The number of `--fq1` options should be same as the number of `--fq2` and `--rg` options.
The order of `--fq1`, `--fq2`, and `--rg` options are important. 

```
$ cd /path/to/working/directory
$ mkdir -p tutorial_02
$ . cwlenv/bin/activate
$ cwltool --singularity \
    --outdir tutorial_02 \
    WGSpipeline/Workflows/germline-gpu.cwl \
    --ref reference_hg38/Homo_sapiens_assembly38.fasta \
    --fq1 wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06HDADXX130110.1\\tPL:ILLUMINA\\tPU:H06HDADXX130110.1\\tLB:H06HDADXX130110.1\\tSM:NA12878" \
    --fq1 wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06HDADXX130110.2\\tPL:ILLUMINA\\tPU:H06HDADXX130110.2\\tLB:H06HDADXX130110.2\\tSM:NA12878" \
    --fq1 wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06JUADXX130110.1\\tPL:ILLUMINA\\tPU:H06JUADXX130110.1\\tLB:H06JUADXX130110.1\\tSM:NA12878" \
    --num_gpus 4 \
    --prefix NA12878 \
    --autosome_interval WGSpipeline/interval_files/autosome.bed \
    --PAR_interval WGSpipeline/interval_files/PAR.bed \
    --chrX_interval WGSpipeline/interval_files/chrX.bed \
    --chrY_interval WGSpipeline/interval_files/chrY.bed
```

Output files will be saved in the directory `/path/to/working/directory/tutorial_02`:
```
/path/to/working/directory/tutorial_02
|--NA12878.PAR.g.vcf.gz
|--NA12878.PAR.g.vcf.gz.tbi
|--NA12878.autosome.g.vcf.gz
|--NA12878.autosome.g.vcf.gz.tbi
|--NA12878.bqsr.recal.table
|--NA12878.chrX_female.g.vcf.gz
|--NA12878.chrX_female.g.vcf.gz.tbi
|--NA12878.chrX_male.g.vcf.gz
|--NA12878.chrX_male.g.vcf.gz.tbi
|--NA12878.chrY.g.vcf.gz
|--NA12878.chrY.g.vcf.gz.tbi
|--NA12878.cram
|--NA12878.cram.crai
```

## Tutorial 3: Run `germline-gpu.cwl` workflow with --knownSites option
Run `germline-gpu.cwl` workflow with --knownSites option

```
$ cd /path/to/working/directory
$ mkdir -p tutorial_03
$ . cwlenv/bin/activate
$ cwltool --singularity \
    --outdir tutorial_03 \
    WGSpipeline/Workflows/germline-gpu.cwl \
    --ref reference_hg38/Homo_sapiens_assembly38.fasta \
    --knownSites reference_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites reference_hg38/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --fq1 wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06HDADXX130110.1\\tPL:ILLUMINA\\tPU:H06HDADXX130110.1\\tLB:H06HDADXX130110.1\\tSM:NA12878" \
    --fq1 wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06HDADXX130110.2\\tPL:ILLUMINA\\tPU:H06HDADXX130110.2\\tLB:H06HDADXX130110.2\\tSM:NA12878" \
    --fq1 wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_1.fastq \
    --fq2 wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_2.fastq \
    --rg "@RG\\tID:NA12878.H06JUADXX130110.1\\tPL:ILLUMINA\\tPU:H06JUADXX130110.1\\tLB:H06JUADXX130110.1\\tSM:NA12878" \
    --num_gpus 4 \
    --prefix NA12878 \
    --autosome_interval WGSpipeline/interval_files/autosome.bed \
    --PAR_interval WGSpipeline/interval_files/PAR.bed \
    --chrX_interval WGSpipeline/interval_files/chrX.bed \
    --chrY_interval WGSpipeline/interval_files/chrY.bed
```

When `--knownSites` option is used, then BQSR will be applied and `_prefix_.bqsr.recla.table` file will not be empty. 

Output files will be saved in the directory `/path/to/working/directory/tutorial_02`:
```
/path/to/working/directory/tutorial_03
|--NA12878.PAR.g.vcf.gz
|--NA12878.PAR.g.vcf.gz.tbi
|--NA12878.autosome.g.vcf.gz
|--NA12878.autosome.g.vcf.gz.tbi
|--NA12878.bqsr.recal.table
|--NA12878.chrX_female.g.vcf.gz
|--NA12878.chrX_female.g.vcf.gz.tbi
|--NA12878.chrX_male.g.vcf.gz
|--NA12878.chrX_male.g.vcf.gz.tbi
|--NA12878.chrY.g.vcf.gz
|--NA12878.chrY.g.vcf.gz.tbi
|--NA12878.cram
|--NA12878.cram.crai
```
