#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//Minimap_BWA
params.vcf_1 = ''
params.vcf_2 = ''
params.output = ''

process BCF_ISEC{
    output:
    stdout

    script:
    """
    if ! [[ -d "${params.output}" ]]; then
    mkdir -p "$baseDir/${params.output}"
    fi
    cd $baseDir
    apptainer exec $baseDir/vip_phase/images/bcftools-1.17.sif /usr/local/bin/bcftools isec -p $params.output $params.vcf_1 $params.vcf_2
    """
}

process BGZIP{
    input:
    stdin

    output:
    stdout

    script:
    """
    cd $baseDir
    apptainer exec $baseDir/images/tabix.sif /opt/samtools/bin/bgzip $params.output/*
    """
}

process BCF_INDEX{
    input:
    stdin

    output:
    stdout

    script:
    """
    cd $baseDir
    apptainer exec $baseDir/vip_phase/images/bcftools-1.17.sif /usr/local/bin/bcftools index $params.output/*
    """
}

workflow{
    BCF_ISEC | BGZIP | BCF_INDEX
}