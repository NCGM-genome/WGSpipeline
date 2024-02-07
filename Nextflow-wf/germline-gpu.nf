#!/usr/bin/env nextflow
params.mode = 'cram'
params.ploidy_1 = 1
params.ploidy_2 = 2

process fq2cram {
    publishDir params.outdir, mode:'copy'

    input:
    val fq1_list
    val fq2_list
    val rg_list
    path ref
    val bwa_options
    val prefix
    val num_gpus
    val mode

    output:
    path "${prefix}.${mode}", emit: cram
    path "${prefix}.bqsr.recal.table"
    path "${prefix}.${mode}.crai"

    script:
    """
    bash /tools/pbrun-fq2bam.sh '$fq1_list' '$fq2_list' '$rg_list' $ref '$bwa_options' $num_gpus $prefix $mode ''
    """
}

process HC_autosome {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path autosome_interval
    val ploidy_2
    val num_gpus
    val HC_autosome_out

    output:
    path "${HC_autosome_out}", emit: gvcf
    path "${HC_autosome_out}.tbi"

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $autosome_interval --ploidy $ploidy_2 --num-gpus $num_gpus --out-variants $HC_autosome_out --gvcf
    """
}

process HC_PAR {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path PAR_interval
    val ploidy_2
    val num_gpus
    val HC_PAR_out

    output:
    path "${HC_PAR_out}", emit: gvcf
    path "${HC_PAR_out}.tbi"

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $PAR_interval --ploidy $ploidy_2 --num-gpus $num_gpus --out-variants $HC_PAR_out --gvcf
    """
}

process HC_chrX_female {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path chrX_interval
    val ploidy_2
    val num_gpus
    val HC_chrX_female_out

    output:
    path "${HC_chrX_female_out}", emit: gvcf
    path "${HC_chrX_female_out}.tbi"

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $chrX_interval --ploidy $ploidy_2 --num-gpus $num_gpus --out-variants $HC_chrX_female_out --gvcf
    """
}

process HC_chrX_male {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path chrX_interval
    val ploidy_1
    val num_gpus
    val HC_chrX_male_out

    output:
    path "${HC_chrX_male_out}", emit: gvcf
    path "${HC_chrX_male_out}.tbi"

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $chrX_interval --ploidy $ploidy_1 --num-gpus $num_gpus --out-variants $HC_chrX_male_out --gvcf
    """
}

process HC_chrY {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path chrY_interval
    val ploidy_1
    val num_gpus
    val HC_chrY_out

    output:
    path "${HC_chrY_out}", emit: gvcf
    path "${HC_chrY_out}.tbi"

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $chrY_interval --ploidy $ploidy_1 --num-gpus $num_gpus --out-variants $HC_chrY_out --gvcf
    """
}

workflow {
    // #### fq2cram ch####
    // 「fastq_reads_1_*」にマッチするparamsの値を抽出してリスト化
    fq1_list = params.findAll { it.key.startsWith('fastq_reads_1') }
        .collect { it.value }
        .sort()
        .join(",")
    fq1_list = Channel.of(fq1_list)
    // 「fastq_reads_2_*」にマッチするparamsの値を抽出してリスト化
    fq2_list = params.findAll { it.key.startsWith('fastq_reads_2') }
        .collect { it.value }
        .sort()
        .join(",")
    fq2_list = Channel.of(fq2_list)
    // rg*にマッチするparamsの値を抽出してリスト化
    rg_list = params.findAll { it.key.startsWith('rg_') }
        .collect { it.value }
        .sort()
        .join(",")
    rg_list = Channel.of(rg_list)

    ref = Channel.fromPath(params.ref)
    bwa_options = Channel.of(params.bwa_options)
    prefix = Channel.of(params.prefix)
    num_gpus = Channel.of(params.num_gpus)
    mode = Channel.of(params.mode)
    // knownSites = Channel.fromPath(params.knownSites)

    // #### HC_autosome ch####
    HC_autosome_out = prefix.map { it + ".autosome.g.vcf.gz" }
    autosome_interval = Channel.fromPath(params.autosome_interval)
    ploidy_1 = Channel.of(params.ploidy_1)
    ploidy_2 = Channel.of(params.ploidy_2)

    // #### HC_PAR ch####
    HC_PAR_out = prefix.map { it + ".PAR.g.vcf.gz" }
    PAR_interval = Channel.fromPath(params.PAR_interval)

    // #### HC_chrX_female ch####
    HC_chrX_female_out = prefix.map { it + ".chrX_female.g.vcf.gz" }
    chrX_interval = Channel.fromPath(params.chrX_interval)

    // #### HC_chrX_male ch####
    HC_chrX_male_out = prefix.map { it + ".chrX_male.g.vcf.gz" }
    chrX_interval = Channel.fromPath(params.chrX_interval)

    // #### HC_chrY ch####
    HC_chrY_out = prefix.map { it + ".chrY.g.vcf.gz" }
    chrY_interval = Channel.fromPath(params.chrY_interval)

    // ################ workflow ################
    // #### fq2cram ####
    fq2cram_out = fq2cram(fq1_list, fq2_list, rg_list, ref, bwa_options, prefix, num_gpus, mode)
    out_cram = fq2cram_out.cram
    // #### HC_autosome ####
    HC_autosome = HC_autosome(ref, out_cram, autosome_interval, ploidy_2, num_gpus, HC_autosome_out)
    out_gvcf_autosome = HC_autosome.gvcf
    // #### HC_PAR####
    HC_PAR = HC_PAR(ref, out_cram, PAR_interval, ploidy_2, num_gpus, HC_PAR_out)
    out_gvcf_PAR = HC_PAR.gvcf
    // #### HC_chrX_female####
    HC_chrX_female = HC_chrX_female(ref, out_cram, chrX_interval, ploidy_2, num_gpus, HC_chrX_female_out)
    out_gvcf_chrX_female = HC_chrX_female.gvcf
    // #### HC_chrX_male####
    HC_chrX_male = HC_chrX_male(ref, out_cram, chrX_interval, ploidy_1, num_gpus, HC_chrX_male_out)
    out_gvcf_chrX_male = HC_chrX_male.gvcf
    // #### HC_chrY####
    HC_chrY = HC_chrY(ref, out_cram, chrY_interval, ploidy_1, num_gpus, HC_chrY_out)
    out_gvcf_chrY = HC_chrY.gvcf
}
