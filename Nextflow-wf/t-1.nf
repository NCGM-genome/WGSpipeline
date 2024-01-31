#!/usr/bin/env nextflow
params.mode = 'cram'
// params.outdir ='result'
params.ploidy = 2

process fq2cram {
    publishDir params.outdir, mode:'copy'

    input:
    path fq1
    path fq2
    val rg
    path ref
    val bwa_options
    val prefix
    val num_gpus
    val mode
    // path knownSites

    output:
    path "${prefix}.${mode}", emit: cram
    path "${prefix}.bqsr.recal.table", emit: recal

    script:
    """
    bash /tools/pbrun-fq2bam.sh $fq1 $fq2 '$rg' $ref '$bwa_options' $num_gpus $prefix $mode ''
    """
}

process HC_autosome {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path autosome_interval
    val ploidy
    val num_gpus
    val HC_autosome_out

    output:
    path "${HC_autosome_out}", emit: gvcf

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $autosome_interval --ploidy $ploidy --num-gpus $num_gpus --out-variants $HC_autosome_out --gvcf
    """
}

process HC_PAR {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path PAR_interval
    val ploidy
    val num_gpus
    val HC_PAR_out

    output:
    path "${HC_PAR_out}", emit: gvcf

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $PAR_interval --ploidy $ploidy --num-gpus $num_gpus --out-variants $HC_PAR_out --gvcf
    """
}

process HC_chrX_female {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path chrX_interval
    val ploidy
    val num_gpus
    val HC_chrX_female_out

    output:
    path "${HC_chrX_female_out}", emit: gvcf

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $chrX_interval --ploidy $ploidy --num-gpus $num_gpus --out-variants $HC_chrX_female_out --gvcf
    """
}

process HC_chrX_male {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path chrX_interval
    val ploidy
    val num_gpus
    val HC_chrX_male_out

    output:
    path "${HC_chrX_male_out}", emit: gvcf

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $chrX_interval --ploidy $ploidy --num-gpus $num_gpus --out-variants $HC_chrX_male_out --gvcf
    """
}

process HC_chrY {
    publishDir params.outdir, mode:'copy'
    label 'haplotypecaller'

    input:
    path ref
    path out_cram
    path chrY_interval
    val ploidy
    val num_gpus
    val HC_chrY_out

    output:
    path "${HC_chrY_out}", emit: gvcf

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $chrY_interval --ploidy $ploidy --num-gpus $num_gpus --out-variants $HC_chrY_out --gvcf
    """
}

workflow {
    // #### fq2cram ch####
    fq1 = Channel.fromPath(params.fq1)
    fq2 = Channel.fromPath(params.fq2)
    rg = Channel.of(params.rg)
    ref = Channel.fromPath(params.ref)
    bwa_options = Channel.of(params.bwa_options)
    prefix = Channel.of(params.prefix)
    num_gpus = Channel.of(params.num_gpus)
    mode = Channel.of(params.mode)
    // knownSites = Channel.fromPath(params.knownSites)
    // #### HC_autosome ch####
    HC_autosome_out = prefix.map { it + ".autosome.g.vcf.gz" }
    autosome_interval = Channel.fromPath(params.autosome_interval)
    ploidy = Channel.of(params.ploidy)
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
    fq2cram_out = fq2cram(fq1, fq2, rg, ref, bwa_options, prefix, num_gpus, mode)
    out_cram = fq2cram_out.cram
    out_recal = fq2cram_out.recal
    // #### HC_autosome ####
    out_gvcf_autosome = HC_autosome(ref, out_cram, autosome_interval, ploidy, num_gpus, HC_autosome_out)
    // #### HC_PAR####
    out_gvcf_PAR = HC_PAR(ref, out_cram, PAR_interval, ploidy, num_gpus, HC_PAR_out)
    // #### HC_chrX_female####
    out_gvcf_chrX_female = HC_chrX_female(ref, out_cram, chrX_interval, ploidy, num_gpus, HC_chrX_female_out)
    // #### HC_chrX_male####
    out_gvcf_chrX_male = HC_chrX_male(ref, out_cram, chrX_interval, ploidy, num_gpus, HC_chrX_male_out)
    // #### HC_chrY####
    out_gvcf_chrY = HC_chrY(ref, out_cram, chrY_interval, ploidy, num_gpus, HC_chrY_out)
}