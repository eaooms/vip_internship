#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
    apptainer exec $baseDir/images/bcftools-1.17.sif /usr/local/bin/bcftools isec -p $params.output $params.vcf_1 $params.vcf_2
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
    ls $params.output/*.vcf | xargs -n1 apptainer exec $baseDir/images/tabix.sif /opt/samtools/bin/bgzip
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
    ls $params.output/*.vcf.gz | xargs -n1 apptainer exec $baseDir/images/bcftools-1.17.sif /usr/local/bin/bcftools index
    """
}

workflow{
    BCF_ISEC | BGZIP | BCF_INDEX
}
