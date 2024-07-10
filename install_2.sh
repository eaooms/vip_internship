#!/bin/bash

# Get the directory of the current script
script_dir=$(dirname "$(readlink -f "$0")")

# Define the base directory where images will be saved
mkdir "$script_dir/images/apptainer_cache"

# Define the base directory where resources will be saved
mkdir "$script_dir/resources/HG003"

export APPTAINER_CACHEDIR="$script_dir/images/apptainer_cache"

# List of Docker images to build
images=(
  "hkubal/clair3:latest"
  "hkubal/clair3-mp"
  "biocontainers/bwa:v0.7.17_cv1"
  "pkrusche/hap.py"
  "mgibio/tabix:1.3.1"

)

# List of GIAB dataset URLs to download
datasets=(
  "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
  "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
  "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
)

# Function to build Docker images using Apptainer
pull_docker_images() {
  for image in "${images[@]}"; do
  # Replace '/' and ':' with '_' to create a valid filename
  sanitized_image=$(echo "$image" | tr '/:' '__')
  # Build the image with Apptainer
  apptainer build "$script_dir/images/$sanitized_image.sif" "docker://$image"
done
}

# Function to download GIAB datasets using wget
download_datasets() {
  for url in "${datasets[@]}"; do
    # Get the filename from the URL
    filename=$(basename "$url")
    
    echo "Downloading dataset: $url"
    echo "Saving as: $script_dir/resources/HG003/$filename"
    # Download the dataset with wget
    wget -P "$script_dir/resources/HG003" "$url"
  done
}

# Execute the functions
pull_docker_images
download_datasets
