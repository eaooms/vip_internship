#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference = "$baseDir/resources/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

//Minimap_BWA
params.ont = ''
params.fastq1 = ''
params.fastq2 = ''
params.output = ''

process BWA{
    output:
    stdout

    script:
    """
    if ! [[ -d "${params.output}" ]]; then
    mkdir -p "$baseDir/${params.output}"
    fi
    cd $baseDir
    apptainer exec $baseDir/images/biocontainers_bwa_v0.7.17_cv1.sif bwa mem -M -t 8 $params.reference $params.fastq1 $params.fastq2| apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools view -bh -o $params.output/Ilmn_unsorted.bam
    """
}

process MINI{
    input:
    stdin

    output:
    stdout

    script:
    """
    cd $baseDir
    $baseDir/images/minimap2/minimap2 -t 8 -a $params.reference $params.ont | apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools view -b -o $params.output/Ont_unsorted.bam
    """
}

process SORTONT{
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools sort -@ 8 -o $params.output/Ont_sorted.bam $params.output/Ont_unsorted.bam
    """
}

process INDEXONT{
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools index -@ 8 $params.output/Ont_sorted.bam
    """
}

process SORTILMN{
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools sort -@ 8 -o $params.output/Ilmn_sorted.bam $params.output/Ilmn_unsorted.bam
    
    """
}

process INDEXILMN{
    input:
    stdin

    output:
    stdout

    script:
    """
    apptainer exec $baseDir/images/samtools-1.17-patch1.sif /usr/local/bin/samtools index -@ 8 $params.output/Ilmn_sorted.bam
    """
}

workflow{
    BWA | MINI | SORTONT | SORTILMN | INDEXONT | INDEXILMN
}
