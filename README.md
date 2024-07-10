# Variant Interpretation Pipeline
VIP is a flexible human variant interpretation pipeline for rare disease using state-of-the-art pathogenicity prediction ([CAPICE](https://github.com/molgenis/capice)) and template-based interactive reporting to facilitate decision support.

![Example Report](docs/img/report_example.png)

## Documentation
VIP documentation is available at this link https://molgenis.github.io/vip/.

> [!TIP]
> Visit <a href="https://vip.molgeniscloud.org/">https://vip.molgeniscloud.org/</a> to analyse your own variants

> [!TIP]
> Preprint now available at <a href="https://doi.org/10.1101/2024.04.11.24305656">medRxiv</a>

## Quick Reference

### Requirements
- [GNU-based Linux](https://en.wikipedia.org/wiki/Linux_distribution#Widely_used_GNU-based_or_GNU-compatible_distributions) (e.g. Ubuntu, [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/about)) with [x86_64](https://en.wikipedia.org/wiki/X86-64) architecture
- Bash ≥ 3.2
- Java ≥ 11
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages) (setuid installation)
- 8GB RAM (an estimate, see also the [documentation](https://molgenis.github.io/vip/get_started/requirements/))
- 220GB disk space

### Installation
```bash
git clone https://github.com/molgenis/vip
bash vip/install.sh
```

### Usage of VIP
```bash
usage: vip -w <arg> -i <arg> -o <arg>
  -w, --workflow <arg>  workflow to execute. allowed values: cram, fastq, gvcf, vcf
  -i, --input    <arg>  path to sample sheet .tsv
  -o, --output   <arg>  output folder
  -c, --config   <arg>  path to additional nextflow .cfg (optional)
  -p, --profile  <arg>  nextflow configuration profile (optional)
  -r, --resume          resume execution using cached results (default: false)
  -s, --stub            quickly prototype workflow logic using process script stubs
  -h, --help            print this message and exit
```

# Proof of concept Clair3-MP tool and phasing data

## How to install VIP and test this branch
### Clone repository 
```bash
git clone https://github.com/eaooms/vip_internship.git
cd vip_internship
```

### Install to download tools and resources
```bash
bash install.sh
bash install_2.sh
```

### Test the Clair3-MP workflow
```bash
usage: clair3mp_workflow.sh [-n <arg> -i <arg> -l <arg> -o <arg>]
  -n, --nanopore             <arg>  path to nanopore file
  -i, --illumina_r1            <arg>  path to fastq_r1 file
  -l, --illumina_r2            <arg>  path to fastq_r2 file
  -b, --bed_file            <arg>  path to BED file
  -o, --output           <arg>  output folder
```

### Test the Phase workflow
```bash
usage: phase_workflow.sh [-n <arg> -i <arg> -o <arg>]
  -n, --nanopore             <arg>  path to nanopore bam file
  -i, --illumina           <arg>  path to illumina bam file
  -o, --output           <arg>  output folder
```

### Reproduce with public data
These scripts can be used to reproduce results from the internship.
Data that wants to be used needs to be in the vip_internship directory

Clair3-MP command
```bash
./nextflow-23.10.0-all Clair3-MP.nf --ont "Sorted and indexed Nanopore file" --ilmn "Illumina Bam file" --ouput "name for output map"
```
Clair3 for both Illumina and Nanopore file command
```bash
./nextflow-23.10.0-all Clair3.nf --ont "Sorted and indexed Nanopore file" --ilmn "Illumina Bam file" --ouput "name for output map"
```
Fastq to BAM for both Illumina and Nanopore files command
```bash
./nextflow-23.10.0-all fastq_to_bam.nf --ont "Nanopore fastq file" --fastq1 "Illumina fastq R1 file" --fastq2 "Illumina fastq R2 file" --ouput "name for output map"
```
Haplotype comparison tool command
```bash
./nextflow-23.10.0-all Happy.nf --vcf "VCF file" --ouput "name for output map"
```
Get all unique variants between VCF's command
```bash
./nextflow-23.10.0-all Unique_in_vcf.nf --vcf_1 "First VCF file" --vcf_2 "Second VCF file" --ouput "name for output map"
```

## Developers
To create the documentation pages:
```
pip install mkdocs mkdocs-mermaid2-plugin
mkdocs serve
```

### License
VIP is an aggregate work of many works, each covered by their own licence(s). For the purposes of determining what you can do with specific works in VIP, this policy should be read together with the licence(s) of the relevant tools. For the avoidance of doubt, where any other licence grants rights, this policy does not modify or reduce those rights under those licences.
