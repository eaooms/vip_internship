#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference = "$baseDir/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

//Workflow_parameters 
params.ont = ''
params.ilmn = ''
params.output = ''

process RUN_CLAIR{
    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/hkubal_clair3_latest.sif /opt/bin/run_clair3.sh \
    --bam_fn=$params.ilmn \
    --ref_fn=$params.reference \
    --model_path /opt/models/ilmn \
    --platform="ilmn" \
    --threads=4 \
    --output=$params.output/Ilmn.vcf
    """
}

process WHATSHAPPHASE{
    input:
    stdin

    output:
    stdout

    script:
    """
    whatshap phase -o $params.output/phased_short_with_long.vcf \
    --ignore-read-groups --reference $params.reference $params.output/Ilmn.vcf $params.ont
    """
}
process BGZIP{
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/mgibio_tabix_1.3.1.sif /opt/samtools/bin/bgzip $params.output/phased_short_with_long.vcf
    """
}

process TABIX{
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/mgibio_tabix_1.3.1.sif /opt/samtools/bin/tabix $params.output/phased_short_with_long.vcf.gz
    """
}


process WHATSHAPHAPLOTAG{
    input:
    stdin

    output:
    stdout

    script:
    """
    whatshap haplotag -o $params.output/phased_haplo.bam \
    --ignore-read-groups --reference $params.reference \
    $params.output/phased_short_with_long.vcf.gz $params.ilmn
    """
}
workflow {
    RUN_CLAIR | WHATSHAPPHASE | BGZIP | TABIX | WHATSHAPHAPLOTAG
}
