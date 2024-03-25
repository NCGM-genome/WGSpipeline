# Workflows
This repository contains following workflows:
- [WGSpipeline Nextflow](#wgspipeline-nextflow)**`(germline-gpu.nf)`:** This workflow calculates alignments (cram) and calls variants (gvcf) from sequence reads (fastq) and reference (fasta) using `parabricks` version 4.0.0. Variants will be output as separate files according to the interval and ploidy.  This workflow is the Nextflow version of [germline-gpu.cwl](../Workflows/germline-gpu.cwl).
- [CNVkit Nextflow](#cnvkit-nextflow)**`(cnvkit.nf)`:** This workflow executes the cnvkit batch command of CNVkit, a Python library and command-line software toolkit for inferring and visualizing copy number variations (CNVs) from DNA sequencing data, using Nextflow.
- [Manta Nextflow](#manta-nextflow)**`(manta.nf)`:** This workflow uses manta:1.6.0 to generate a configuration file (configManta.py) and execute the workflow (runWorkflow.py), enabling the fast and accurate detection of structural variations (SVs) and indels, using Nextflow.
- [Sentieon joint calling Nextflow](#sentieon-joint-calling-nextflow)**`(sentieon-jc.nf)`:** This workflow executes the sentioen command with the gvcf file as input, using Nextflow.

# WGSpipeline Nextflow
## Installation requirements
- [Hardware requirements to run parabricks](https://docs.nvidia.com/clara/parabricks/4.0.0/GettingStarted.html#hardware-requirements)
- [Software requirements to run parabricks](https://docs.nvidia.com/clara/parabricks/4.0.0/GettingStarted.html#software-requirements)
- SingularityCE 4.0.0+
- openjdk 11.0.20.1+
- Nextflow 23.10.1+
## Install `Nextflow`
You can install Nextflow by executing the following commands:
- Check prerequisites：Java 11 or later is required
  ```
  $ java -version
  ```
- Set up：Dead easy to install
  ```
  $ curl -s https://get.nextflow.io | bash
  ```
- Launch：Try a simple demo
  ```
  $ ./nextflow run hello
  ```

- Set path
  ```
  export PATH=$PATH:/path/to/nextflow/directory
  ```

Please confirm the nextflow version:
```
  $ nextflow info
  Version: 23.10.1 build 5891
  Created: 12-01-2024 22:01 UTC (13-01-2024 07:01 JDT)
  System: Linux 5.15.0-88-generic
  Runtime: Groovy 3.0.19 on OpenJDK 64-Bit Server VM 11.0.20.1+1-post-Ubuntu-0ubuntu122.04
  Encoding: UTF-8 (UTF-8)
```
Install guide：[Getting started](https://www.nextflow.io/)
## Clone this repository
Clone this repository by executing the following commands:
```
$ cd /path/to/working/directory/
$ git clone https://github.com/NCGM-genome/WGSpipeline.git
```
## Preparation of `germline-gpu.nf` workflow
- Change Nextflow-wf directory
  ```
  cd WGSpipeline/Nextflow-wf
  ```
- Creation of input file
  ```
  touch example.config
  ```
## Usage of `germline-gpu.nf` workflow
```
nextflow run germline-gpu.nf -c example.config
```
## Input file
- example.config
```groovy
singularity {
    enabled = true
}

process {
    withName: fq2cram {
        containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
        container = 'docker://hacchy/pbrun-fq2bam:4.0.0-1_v20230412'
        queue = '<slurm partition name>'
        executor = 'slurm'
        memory = '<Maximum memory value> GB'
    }
}

process { 
    withLabel: haplotypecaller { 
        containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
        container = 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1'
        queue = '<slurm partition name>'
        executor = 'slurm'
        memory = '<Maximum memory value> GB'
    }
}

params {
    // Output directory.
    outdir = 'path/to/dir'
    // Path to FASTQ file 1. FASTQ files defined by "fastq_reads_1_*". Please write the full path.
    fastq_reads_1_* = 'path/to/file'
    // Path to FASTQ file 2. FASTQ files defined by "fastq_reads_2_*". Please write the full path.
    fastq_reads_2_* = 'path/to/file'
    // Read group string. Read group defined by "rg_*"
    rg_* = 'STRING'
    // Path to the reference file.
    ref = 'path/to/file'
    // Pass supported bwa mem options as one string.
    bwa_options = '-T 0 -Y'
    // Output file prefix.
    prefix = 'STRING'
    // Number of GPUs to use for a run (should be ≥1).
    num_gpus = INT
    // Path to interval BED file for autosome regions.
    autosome_interval = 'path/to/file'
    // Path to interval BED file for PAR regions.
    PAR_interval = 'path/to/file'
    // Path to interval BED file for chrX regions.
    chrX_interval = 'path/to/file'
    // Path to interval BED file for chrY regions.
    chrY_interval = 'path/to/file'
    // Path to a known indels file. knownSites files defined by "knownSites_*".
    knownSites_* = 'path/to/file'
}
```
- This config file assumes the following execution conditions
  - singularity as container runtime
  - slurm as executor
  - GPU Nodes
- Memo
  - The **`slurm partition name`:**  can be checked with **`sinfo -l`:** 
## Download reference and resource files
Download reference and resource files from the URLs listed in [reference_hg38.download_links.txt](./download_links/reference_hg38.download_links.txt) by executing the following commnds:
```
$ OUTDIR=reference_hg38 ; mkdir -p $OUTDIR ; for url in `cat ../download_links/reference_hg38.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```
## Download a dataset for tutorial
Download a dataset for tutorial from the URLs listed in [wgs_fastq_NA12878_20k.download_links.txt](./download_links/wgs_fastq_NA12878_20k.download_links.txt) by executing the following commnds:
```
$ OUTDIR=wgs_fastq ; mkdir -p $OUTDIR ; for url in `cat ../download_links/wgs_fastq_NA12878_20k.download_links.txt` ; do echo $url ; file=`basename $url` ; if [ ! -f ${OUTDIR}/$file ] ; then wget $url -O ${OUTDIR}/$file ; fi ; done
```
## Tutorial 1: Run  `germline-gpu.nf` workflow with a pair of FASTQ files.
Run `germline-gpu.nf` workflow with a pair of FASTQ files. 
- execution
  ```
  nextflow run germline-gpu.nf -c tutorial_01.config
  ```
- tutorial_01.config
  ```groovy
  singularity {
      enabled = true
  }

  process {
      withName: fq2cram {
          containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
          container = 'docker://hacchy/pbrun-fq2bam:4.0.0-1_v20230412'
          queue = '<slurm partition name>'
          executor = 'slurm'
          memory = '<Maximum memory value> GB'
      }
  }

  process { 
      withLabel: haplotypecaller { 
          containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
          container = 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1'
          queue = '<slurm partition name>'
          executor = 'slurm'
          memory = '<Maximum memory value> GB'
      }
  }

  params {
      outdir = 'tutorial_01'
      fastq_reads_1_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_2_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq'
      rg_1 = '@RG\\tID:NA12878.H06HDADXX130110.1\\tPL:ILLUMINA\\tPU:H06HDADXX130110.1\\tLB:H06HDADXX130110.1\\tSM:NA12878'
      ref = 'reference_hg38/Homo_sapiens_assembly38.fasta'
      bwa_options = '-T 0 -Y'
      prefix = 'NA12878.H06HDADXX130110.1'
      num_gpus = 4
      autosome_interval = '../interval_files/autosome.bed'
      PAR_interval = '../interval_files/PAR.bed'
      chrX_interval = '../interval_files/chrX.bed'
      chrY_interval = '../interval_files/chrY.bed'
  }
  ```
When `knownSites` option is not used, then empty BQSR table file `<prefix>.bqsr.recla.table` will be created. 

Output files will be saved in the directory `/path/to/working/directory/WGSpipeline/Nextflow-wf/tutorial_01`:
```
/path/to/working/directory/WGSpipeline/Nextflow-wf/tutorial_01
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
## Tutorial 2: Run `germline-gpu.nf` workflow with multiple pairs of FASTQ files.
Run `germline-gpu.nf` workflow with multiple pairs of FASTQ files.

`fastq_reads_1_*`, `fastq_reads_2_*`, and `rg_*` params can be repeated multiple times. 
The number of `fastq_reads_1_*` params should be same as the number of `fastq_reads_2_*` and `rg_*` params.
- execution
  ```
  nextflow run germline-gpu.nf -c tutorial_02.config
  ```
- tutorial_02.config
  ```groovy
  singularity {
      enabled = true
  }

  process {
      withName: fq2cram {
          containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
          container = 'docker://hacchy/pbrun-fq2bam:4.0.0-1_v20230412'
          queue = '<slurm partition name>'
          executor = 'slurm'
          memory = '<Maximum memory value> GB'
      }
  }

  process { 
      withLabel: haplotypecaller { 
          containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
          container = 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1'
          queue = '<slurm partition name>'
          executor = 'slurm'
          memory = '<Maximum memory value> GB'
      }
  }

  params {
      outdir = 'tutorial_02'
      fastq_reads_1_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_1_2 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_1_3 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_2_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq'
      fastq_reads_2_2 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_2.fastq'
      fastq_reads_2_3 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_2.fastq'
      rg_1 = '@RG\\tID:NA12878.H06HDADXX130110.1\\tPL:ILLUMINA\\tPU:H06HDADXX130110.1\\tLB:H06HDADXX130110.1\\tSM:NA12878'
      rg_2 = '@RG\\tID:NA12878.H06HDADXX130110.2\\tPL:ILLUMINA\\tPU:H06HDADXX130110.2\\tLB:H06HDADXX130110.2\\tSM:NA12878'
      rg_3 = '@RG\\tID:NA12878.H06JUADXX130110.1\\tPL:ILLUMINA\\tPU:H06JUADXX130110.1\\tLB:H06JUADXX130110.1\\tSM:NA12878'
      ref = 'reference_hg38/Homo_sapiens_assembly38.fasta'
      bwa_options = '-T 0 -Y'
      prefix = 'NA12878'
      num_gpus = 4
      autosome_interval = '../interval_files/autosome.bed'
      PAR_interval = '../interval_files/PAR.bed'
      chrX_interval = '../interval_files/chrX.bed'
      chrY_interval = '../interval_files/chrY.bed'
  }
  ```

Output files will be saved in the directory `/path/to/working/directory/WGSpipeline/Nextflow-wf/tutorial_02`:
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
## Tutorial 3: Run `germline-gpu.nf` workflow with --knownSites option
Run `germline-gpu.nf` workflow with knownSites params.

- execution
  ```
  nextflow run germline-gpu.nf -c tutorial_03.config
  ```
- tutorial_03.config
  ```groovy
  singularity {
      enabled = true
  }

  process {
      withName: fq2cram {
          containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
          container = 'docker://hacchy/pbrun-fq2bam:4.0.0-1_v20230412'
          queue = '<slurm partition name>'
          executor = 'slurm'
          memory = '<Maximum memory value> GB'
      }
  }

  process { 
      withLabel: haplotypecaller { 
          containerOptions = '--nv --bind /path/to/cuda:/usr/local/cuda'
          container = 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1'
          queue = '<slurm partition name>'
          executor = 'slurm'
          memory = '<Maximum memory value> GB'
      }
  }

  params {
      outdir = 'tutorial_03'
      fastq_reads_1_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_1_2 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_1_3 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_1.fastq'
      fastq_reads_2_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq'
      fastq_reads_2_2 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06HDADXX130110.2.ATCACGAT.20k_reads_2.fastq'
      fastq_reads_2_3 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/wgs_fastq/H06JUADXX130110.1.ATCACGAT.20k_reads_2.fastq'
      rg_1 = '@RG\\tID:NA12878.H06HDADXX130110.1\\tPL:ILLUMINA\\tPU:H06HDADXX130110.1\\tLB:H06HDADXX130110.1\\tSM:NA12878'
      rg_2 = '@RG\\tID:NA12878.H06HDADXX130110.2\\tPL:ILLUMINA\\tPU:H06HDADXX130110.2\\tLB:H06HDADXX130110.2\\tSM:NA12878'
      rg_3 = '@RG\\tID:NA12878.H06JUADXX130110.1\\tPL:ILLUMINA\\tPU:H06JUADXX130110.1\\tLB:H06JUADXX130110.1\\tSM:NA12878'
      ref = 'reference_hg38/Homo_sapiens_assembly38.fasta'
      bwa_options = '-T 0 -Y'
      prefix = 'NA12878'
      num_gpus = 4
      autosome_interval = '../interval_files/autosome.bed'
      PAR_interval = '../interval_files/PAR.bed'
      chrX_interval = '../interval_files/chrX.bed'
      chrY_interval = '../interval_files/chrY.bed'
      knownSites_1 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/reference_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
      knownSites_2 = '/path/to/working/directory/WGSpipeline/Nextflow-wf/reference_hg38/Homo_sapiens_assembly38.known_indels.vcf.gz'
  }
  ```

  When `knownSites` params is used, then non-empty BQSR table file `<prefix>.bqsr.recla.table` will be created. 
  Note that BQSR is not applied to the output cram file. 

  Output files will be saved in the directory `/path/to/working/directory/WGSpipeline/Nextflow-wf/tutorial_03`:
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

# CNVkit Nextflow
## Installation requirements
- SingularityCE 4.0.0+
- openjdk 11.0.20.1+
- Nextflow 23.10.1+
## [Install `Nextflow`](#install-nextflow)
- See previous section.
## [Clone this repository](#clone-this-repository)
- See previous section.
## Preparation of `cnvkit.nf` workflow
- Change Nextflow-wf directory
  ```
  cd WGSpipeline/Nextflow-wf
  ```
- Creation of input file
  ```
  touch example.config
  ```
## Usage of `cnvkit.nf` workflow
```
nextflow run cnvkit.nf -c example.config
```
## Input file
- example.config
```groovy
singularity {
    enabled = true
}

process {
    withLabel: samtools {
        container = 'docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'
        clusterOptions = "--mem-per-cpu 8000M"
        executor = 'slurm'
        queue = '<slurm partition name>'
        cpus = '<cpu core number>'
    }
}

process {
    withName: cnvkit_batch {
        container = 'docker://quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0'
        executor = 'slurm'
        queue = '<slurm partition name>'
        clusterOptions = "--mem-per-cpu 8000M"
        cpus = '<cpu core number>'
    }
}

params {
    // Output directory.
    outdir = 'path/to/dir'
    // Path to cram list. List file containing paths to CRAM files on each line
    sample_list = 'path/to/file'
    // Path to reference file(.fasta)
    ref = 'path/to/file'
    // Path to reference file index(.fasta.fai)
    ref_idx = 'path/to/file'
    // Path to .reference.cnn
    cnn =  'path/to/file'
    // Number of samtools_thread. samtools_thread < cpus (withLabel: samtools)
    samtools_thread = '<thread number>'
    // Number of ncore. ncore < cpus  (withName: cnvkit_batch)
    ncore = '<cnvkit batch cpu core number>'
}
```
- This config file assumes the following execution conditions
  - singularity as container runtime
  - slurm as executor
  - CPU Nodes
- Memo
  - The **`slurm partition name`:**  can be checked with **`sinfo -l`:** 

# Manta Nextflow
## Installation requirements
- SingularityCE 4.0.0+
- openjdk 11.0.20.1+
- Nextflow 23.10.1+
## [Install `Nextflow`](#install-nextflow)
- See previous section.
## [Clone this repository](#clone-this-repository)
- See previous section.
## Preparation of `manta.nf` workflow
- Change Nextflow-wf directory
  ```
  cd WGSpipeline/Nextflow-wf
  ```
- Creation of input file
  ```
  touch example.config
  ```
## Usage of `manta.nf` workflow
```
nextflow run manta.nf -c example.config
```
## Input file
- example.config
```groovy
singularity {
    enabled = true
}

process {
    withLabel: samtools {
        container = 'docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'
        clusterOptions = "--mem-per-cpu 8000M"
        executor = 'slurm'
        queue = '<slurm partition name>'
        cpus = '<cpu core number>'
    }
}

process {
    withName: configManta {
        container = 'docker://fredhutch/manta:1.6.0'
        executor = 'slurm'
        queue = '<slurm partition name>'
        clusterOptions = "--mem-per-cpu 8000M"
        cpus = '<cpu core number>'
    }
}

params {
    // Output directory.
    outdir = 'path/to/dir'
    // Path to cram list. List file containing paths to CRAM files on each line
    sample_list = 'path/to/file'
    // Path to reference file(.fasta)
    ref = 'path/to/file'
    // Path to reference file index(.fasta.fai)
    ref_idx = 'path/to/file'
    // Path to region　file(.bed.gz)
    region = '/lustre8/home/yamaken-gaj-pg/Projects/manta/manta/region.bed.gz'
    // Path to region　file index(.bed.gz.tbi)
    region_idx = '/lustre8/home/yamaken-gaj-pg/Projects/manta/manta/region.bed.gz.tbi'
    // Number of samtools_thread. samtools_thread < cpus (withLabel: samtools)
    samtools_thread = '<thread number>'
    // Number of ncore. ncore < cpus  (withName: configManta)
    ncore = '<runWorkflow.py cpu core number>'
}
```
- This config file assumes the following execution conditions
  - singularity as container runtime
  - slurm as executor
  - CPU Nodes
- Memo
  - The **`slurm partition name`:**  can be checked with **`sinfo -l`:** 

# Sentieon joint calling Nextflow
## Installation requirements
- openjdk 11.0.20.1+
- Nextflow 23.10.1+
- jemalloc 5.3.0-147-ge4817c8d89a2a413e835c4adeab5c5c4412f9235
- sentieon sentieon-genomics-202308
## Install `jemalloc`
- Install jemalloc to use with Sentieon
  ```
  $ git clone https://github.com/jemalloc/jemalloc.git
  $ cd jemalloc
  $ ./autogen.sh
  $ ./configure --prefix=$HOME/workdir/jemalloc_install
  $ make
  $ make install
  $ export PATH=$HOME/workdir/jemalloc_install/bin:$PATH
  ```
## Install `sentieon`
- Get `sentieon-genomics-202308.tar.gz`, unzip, add sentieon to your path
  ```
  $ tar zxvf sentieon-genomics-202308.tar.gz
  $ export PATH=path/to/workdir/sentieon-genomics-202308/bin:$PATH
  ```
- Set up Sentieon license server
  ```
  $ export SENTIEON_LICENSE = <your license server>
  ```
## [Install `Nextflow`](#install-nextflow)
- See previous section.
## [Clone this repository](#clone-this-repository)
- See previous section.
## Preparation of `sentieon-jc.nf` workflow
- Change Nextflow-wf directory
  ```
  cd WGSpipeline/Nextflow-wf
  ```
- Creation of input file
  ```
  touch example.config
  ```
## Usage of `sentieon-jc.nf` workflow
```
nextflow run sentieon-jc.nf -c example.config
```
## Input file
- example.config
```groovy
singularity {
    enabled = true
}

// delete work directory files
cleanup = true

params {
    // Output directory.
    outdir = 'path/to/dir'
    // Prefix name
    prefix = 'PREFIX'
    // Path to gvcf list (gvcfs.txt). List file containing paths to gvcf files on each line
    gvcfs_list = 'path/to/file'
    // Shard options List (shard.txt). List file containing options on each line
    shard = 'path/to/file'
    // Path to reference file(.fasta)
    ref = 'path/to/file'
    // Path to reference file index(.fasta.fai)
    ref_idx = 'path/to/file'
    // Number of ncore.
    ncore = 64
    // Directory of Sentieon commands
    sentieon = 'path/to/dir'
    // Parallel processing information for each chromosome (shard_chr.txt)
    shard_chr = 'path/to/file'
    // Directory of input VCF and its index files（dbsnp, hapmap, omni, i1000g, mills, axiom）
    bundle = 'path/to/dir'
    // VCF dbsnp
    dbsnp = 'dbsnp_151.hg38.vcf.gz'
    // VCF hapmap
    hapmap = 'hapmap_3.3.hg38.vcf.gz'
    // VCF omni
    omni = '1000G_omni2.5.hg38.vcf.gz'
    // VCF i1000g
    i1000g = '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
    // VCF mills
    mills = 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
    // VCF axiom
    axiom = 'Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz'
}

process {
    withLabel: bamtols {
        container = '/path/to/bcftools:1.16--hfe4b78e_1'
        containerOptions = '--bind /path/to/home:/path/to/home'
    }
}

process {
    withLabel: tabix {
        container = '/path/to/bcftools:1.10.2--hd2cd319_0'
        containerOptions = '--bind /path/to/home:/path/to/home'
    }
}
```
- This config file assumes the following execution conditions
  - singularity as container runtime
  - CPU Nodes
