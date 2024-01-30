#!/usr/bin/env nextflow
params.mode = 'cram'
// params.outdir ='result'
params.ploidy = 2

process fq2cram {
    publishDir params.outdir, mode:'copy'
    label 'gpu'

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
    label 'gpu'

    input:
    path ref
    path out_cram
    path interval_file
    val ploidy
    val num_gpus
    val HC_autosome_out

    output:
    path "${HC_autosome_out}", emit: gvcf

    script:
    """
    pbrun haplotypecaller --ref $ref --in-bam $out_cram --interval-file $interval_file --ploidy $ploidy --num-gpus $num_gpus --out-variants $HC_autosome_out --gvcf
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
    interval_file = Channel.fromPath(params.interval_file)
    ploidy = Channel.of(params.ploidy)

    // ######## workflow ########
    // #### fq2cram ####
    fq2cram_out = fq2cram(fq1, fq2, rg, ref, bwa_options, prefix, num_gpus, mode)
    out_cram = fq2cram_out.cram
    out_recal = fq2cram_out.recal
    // #### HC_autosome ####
    out_gvcf_autosome = HC_autosome(ref, out_cram, interval_file, ploidy, num_gpus, HC_autosome_out)
}