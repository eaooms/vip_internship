#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference = "$baseDir/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

//Minimap_BWA
params.ont = ''
params.ilmn = ''
params.output = ''

process RUN_CLAIR_ONT{
    output:
    stdout

    script:
    """
    if ! [[ -d "${params.output}" ]]; then
    mkdir -p "$baseDir/${params.output}"
    fi
    cd $baseDir
    apptainer exec $baseDir/images/clairLatest.sif /opt/bin/run_clair3.sh \
    --bam_fn=$params.ont \
    --ref_fn=$params.reference \
    --model_path /opt/models/r941_prom_sup_g5014 \
    --platform="ont" \
    --threads=4 \
    --output=$baseDir/$params.output/Nanopore.vcf
    """
}

process RUN_CLAIR_ILMN{
    input:
    stdin

    output:
    stdout

    script:
    """
    cd $baseDir
    apptainer exec $baseDir/images/clairLatest.sif /opt/bin/run_clair3.sh \
    --bam_fn=$params.ilmn \
    --ref_fn=$params.reference \
    --model_path /opt/models/ilmn \
    --platform="ilmn" \
    --threads=4 \
    --output=$params.output/Ilmn.vcf
    """
}

workflow{
    RUN_CLAIR_ONT | RUN_CLAIR_ILMN
}