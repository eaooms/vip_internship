#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "$0")")
SCRIPT_NAME="$(basename "$0")"

# SCRIPT_DIR is incorrect when vip.sh is submitted as a Slurm job that is submitted as part of another Slurm job
WORK_DIR="${WORK_DIR:-"${SCRIPT_DIR}"}"

usage() {
  echo -e "usage: ${SCRIPT_NAME} [-n <arg> -i <arg> -l <arg> -o <arg>]
  -n, --nanopore             <arg>  path to nanopore file
  -i, --illumina_r1            <arg>  path to fastq_r1 file
  -l, --illumina_r2            <arg>  path to fastq_r2 file
  -b, --bed_file            <arg>  path to BED file
  -o, --output           <arg>  output folder"
}

validate() {
  local -r nanopore="${1}"
  local -r illumina_r1="${2}"
  local -r illumina_r2="${3}"
  local -r output="${4}"
  local -r bed_file="${5}"

  
  if [[ -z "${nanopore}" ]]; then
    >&2 echo -e "error: missing required -n / --input"
    usage
    exit 2
  fi
  if [[ ! -f "${nanopore}" ]] ; then
    >&2 echo -e "error: input '${nanopore}' does not exist"
    exit 2
  fi

  if [[ -z "${illumina_r1}" ]]; then
    >&2 echo -e "error: missing required -i / --input"
    usage
    exit 2
  fi
  if [[ ! -f "${illumina_r1}" ]] ; then
    >&2 echo -e "error: input '${illumina_r1}' does not exist"
    exit 2
  fi

  if [[ -z "${illumina_r2}" ]]; then
    >&2 echo -e "error: missing required -l / --input"
    usage
    exit 2
  fi
  if [[ ! -f "${illumina_r2}" ]] ; then
    >&2 echo -e "error: input '${illumina_r2}' does not exist"
    exit 2
  fi

  if [[ -z "${bed_file}" ]]; then
    >&2 echo -e "error: missing required -b / --input"
    usage
    exit 2
  fi
  if [[ ! -f "${bed_file}" ]] ; then
    >&2 echo -e "error: input '${bed_file}' does not exist"
    exit 2
  fi

  if [[ -z "${output}" ]]; then
    >&2 echo -e "error: missing required -o / --output"
    usage
    exit 2
  fi

  # detect java, try to load module with name 'java' or 'Java' otherwise
  if ! command -v java &> /dev/null; then
    if command -v module &> /dev/null; then
      if module is_avail java; then
        module load java
      elif module is_avail Java; then
        module load Java
      else
        >&2 echo -e "error: missing required 'java'. could not find a module with name 'java' or 'Java' to load"
        exit 2
      fi
    else
      >&2 echo -e "error: missing required 'java'"
      exit 2
    fi
  fi
}

execute_workflow(){
  local -r paramNanopore="$(realpath "${1}")"
  local -r paramIllumina_r1="$(realpath "${2}")"
  local -r paramIllumina_r2="$(realpath "${3}")"
  local -r paramOutput="$(realpath "${4}")"
  local -r ID_ont="${1}"
  local -r ID_illumina="${2}"
  local -r paramProfile="${5}"
  local -r bedFile="${6}"


  $WORK_DIR/nextflow-23.10.0-all clair3mp_workflow.nf -c nextflow.config -profile $paramProfile --ont $paramNanopore --fastq1 $paramIllumina_r1 --fastq2 $paramIllumina_r2 --output $paramOutput --nanopore_id $ID_ont --ilmn_id $ID_illumina --bed_file $bedFile
}

main() {
  local nanopore=""
  local illumina_r1=""
  local illumina_r2=""
  local bed_file=""
  local output=""
  local profile=""
  if command -v sbatch &> /dev/null; then
    profile="slurm"
  else
    profile="local"
  fi
  local args=$(getopt -o n:i:l:b:o: --long nanopore:,illumina_r1:,illumina_r2:,bed_file:,output:,help: -- "$@")
  eval set -- "${args}"

  while :; do
    case "$1" in
    -n | --nanopore)
      nanopore="$2"
      shift 2 ;;
    -i | --illumina_r1 )
      illumina_r1="$2"
      shift 2 ;;
    -l | --illumina_r2)
      illumina_r2="$2"
      shift 2 ;;
    -b | --bed_file)
      bed_file="$2"
      shift 2 ;;
    -o | --output)
      output="$2"
      shift 2 ;;
    --)
      shift
      break ;;
    *)
      echo "Invalid option : $1"
      exit 1 ;;
    esac
  done

  validate "${nanopore}" "${illumina_r1}" "${illumina_r2}" "${output}" "${bed_file}" 

  if ! [[ -d "${output}" ]]; then
    mkdir -p "${output}"
  fi

  execute_workflow "${nanopore}" "${illumina_r1}" "${illumina_r2}" "${output}" "${profile}" "${bed_file}" 
}

main "${@}"
