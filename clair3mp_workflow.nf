#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters with default values
params.reference = "$baseDir/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
params.scriptpy = "$baseDir/images/tsvoutput.py"

// Inputs
params.ont = ''
params.fastq1 = ''
params.fastq2 = ''
params.output = ''
params.nanopore_id = ''
params.ilmn_id = ''
params.bed_file = ''

// Extract sample IDs from input filenames
params.nanoporeID = params.nanopore_id.split('\\.')[0]
params.ilmnID = params.ilmn_id.split('\\.')[0]

// BWA process for aligning Illumina reads
process BWA {
    output:
    stdout

    script:
    """
    mkdir -p $params.output/intermediates
    apptainer exec $baseDir/images/biocontainers_bwa_v0.7.17_cv1.sif bwa mem -M -t 8 $params.reference $params.fastq1 $params.fastq2 | apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools view -bh -o $params.output/intermediates/Ilmn_unsorted.bam
    """
}

// Minimap2 process for aligning Nanopore reads
process MINI {
    input:
    stdin

    output:
    stdout

    script:
    """
    $baseDir/images/minimap2/minimap2 -t 8 -a $params.reference $params.ont | apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools view -b -o $params.output/intermediates/Ont_unsorted.bam
    """
}

// Process for sorting Nanopore BAM file
process SORTONT {
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools sort -@ 8 -o $params.output/intermediates/Ont_sorted.bam $params.output/intermediates/Ont_unsorted.bam
    """
}

// Process for indexing sorted Nanopore BAM file
process INDEXONT {
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools index -@ 8 $params.output/intermediates/Ont_sorted.bam
    """
}

// Process for sorting Illumina BAM file
process SORTILMN {
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools sort -@ 8 -o $params.output/intermediates/Ilmn_sorted.bam $params.output/intermediates/Ilmn_unsorted.bam
    """
}

// Process for indexing sorted Illumina BAM file
process INDEXILMN {
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools index -@ 8 $params.output/intermediates/Ilmn_sorted.bam
    """
}

// Process for running Clair3 multi-platform variant caller
process CLAIR {
    publishDir "$params.output"

    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec --bind /groups $baseDir/images/clair3mp.sif /opt/bin/run_clair3_mp.sh \
    --bam_fn_c=$params.output/intermediates/Ont_sorted.bam \
    --bam_fn_p1=$params.output/intermediates/Ilmn_sorted.bam \
    --bam_fn_c_platform=ont \
    --bam_fn_p1_platform=ilmn \
    --bed_fn=$params.bed_file \
    --threads=8 \
    --output=$params.output \
    --ref_fn=$params.reference \
    --model_path_clair3_c /opt/bin/models/clair3_models/ont_guppy5 \
    --model_path_clair3_p1 /opt/bin/models/clair3_models/ilmn \
    --model_path_clair3_mp /opt/bin/models/clair3_mp_models/ont_ilmn \
    --sample_name_c=$params.nanoporeID \
    --sample_name_p1=$params.ilmnID
    """
}

// Process for creating TSV output file from VCF
process MAKETSV {
    publishDir "$params.output", mode: 'move'

    input:
    stdin

    output:
    path 'vcf_for_vip.tsv'

    script:
    """
    #!/usr/bin/env python
    import subprocess

    subprocess.run(['python', '$params.scriptpy', 'vcf_for_vip.tsv', '$params.nanoporeID', '$params.output/${params.nanoporeID}.gz'])
    """
}

// Workflow definition, specifying the execution order of processes
workflow {
    BWA | MINI | SORTONT | SORTILMN | INDEXONT | INDEXILMN | CLAIR | MAKETSV
}
