#!/usr/bin/env nextflow

process PbrunFq2CramMultiRGs {
    input:
    path fq1
    path fq2
    val rg
    path ref

    output:
    path "*.cram"

    script:
    """
    singularity exec --nv docker://hacchy/pbrun-fq2bam:4.0.0-1_v20230412 bash /tools/pbrun-fq2bam.sh \
        --fq1 $fq1 \
        --fq2 $fq2 \
        --rg $rg \
        --ref $ref
    """
}

process PbrunHaplotypeCallerFromCram {
    input:
    path cram_file

    output:
    path "*.vcf"

    script:
    """
    singularity exec --nv docker://nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1 pbrun haplotypecaller \
        --input-cram $cram_file \
        --reference $params.ref
    """
}

workflow {
    fq1 = Channel.fromPath(params.fq1)
    fq2 = Channel.fromPath(params.fq2)
    rg = Channel.value(params.rg)
    ref = Channel.fromPath(params.ref)

    cram_files = PbrunFq2CramMultiRGs(fq1, fq2, rg, ref)
    PbrunHaplotypeCallerFromCram(cram_files)
}
