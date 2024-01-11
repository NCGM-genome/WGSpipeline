#!/usr/bin/env nextflow

process Fq2Cram {
    input:
    path fq1
    path fq2
    val rg
    path ref
    path knownSites
    val bwa_options
    val prefix
    val num_gpus

    output:
    path "${prefix}.cram"
    path "${prefix}.cram.crai"
    path "${prefix}.bqsr.recal.table"

    script:
    """
    bash /tools/pbrun-fq2bam.sh \
        ${fq1.join(',')} \
        ${fq2.join(',')} \
        ${rg.join(',')} \
        $ref \
        $bwa_options \
        $num_gpus \
        $prefix \
        $mode
        ${knownSites_ch.collect().join(',')}
    """
}

workflow {
    // 必要な入力ファイルやパラメータを定義
    fq1 = Channel.fromPath(params.fq1)
    fq2 = Channel.fromPath(params.fq2)
    rg = params.rg
    ref = Channel.fromPath(params.ref)
    bwa_options = params.bwa_options ?: "-T 0 -Y"
    num_gpus = params.num_gpus
    prefix = params.prefix
    mode = params.mode
    knownSites_ch = params.knownSites ? Channel.fromPath(params.knownSites) : Channel.empty()

    // Fq2Cram プロセスを呼び出し
    Fq2Cram(fq1, fq2, rg, ref, bwa_options, num_gpus, prefix, mode, knownSites_ch)
}
