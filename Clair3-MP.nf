#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference = "$baseDir/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
params.ont = ''
params.ilmn = ''
params.output = ''

process CLAIR{
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
    apptainer exec --bind /groups $baseDir/images/hkubal_clair3-mp.sif /opt/bin/run_clair3_mp.sh\
    --bam_fn_c=$params.ont\
    --bam_fn_p1=$params.ilmn\
    --bam_fn_c_platform=ont\
    --bam_fn_p1_platform=ilmn\
    --threads=8\
    --output=$baseDir/$params.output\
    --ref_fn=$params.reference\
    --model_path_clair3_c /opt/bin/models/clair3_models/ont_guppy5\
    --model_path_clair3_p1 /opt/bin/models/clair3_models/ilmn\
    --model_path_clair3_mp /opt/bin/models/clair3_mp_models/ont_ilmn\
    --sample_name_c=Ont\
    --sample_name_p1=Ilmn
    """
}

workflow{
    CLAIR(params.reference)
}
