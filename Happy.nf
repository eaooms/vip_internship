#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference = "$baseDir/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
params.truth_vcf = "$baseDir/resources/HG003/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
params.truth_bed = "$baseDir/resources/HG003/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

//Minimap_BWA
params.vcf = ''
params.output = ''

process HAPPY{
    input:
    stdin

    output:
    stdout

    script:
    """
    if ! [[ -d "${params.output}" ]]; then
    mkdir -p "$baseDir/${params.output}"
    fi
    cd $baseDir
    apptainer exec $baseDir/images/pkrusche_hap.py.sif /opt/hap.py/bin/hap.py $params.truth_vcf $params.vcf -f $params.truth_bed -r $params.reference --threads 4 -o $params.output/outputPrefix 
    """
}

workflow{
    HAPPY(params.reference)
}
