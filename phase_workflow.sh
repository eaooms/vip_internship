#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "$0")")
SCRIPT_NAME="$(basename "$0")"

# SCRIPT_DIR is incorrect when vip.sh is submitted as a Slurm job that is submitted as part of another Slurm job
WORK_DIR="${WORK_DIR:-"${SCRIPT_DIR}"}"

usage() {
  echo -e "usage: ${SCRIPT_NAME} [-n <arg> -i <arg> -o <arg>]
  -n, --nanopore             <arg>  path to nanopore bam file
  -i, --illumina           <arg>  path to illumina bam file
  -o, --output           <arg>  output folder"
}

validate() {
  local -r nanopore="${1}"
  local -r illumina="${2}"
  local -r output="${3}"

  if [[ -z "${nanopore}" ]]; then
    >&2 echo -e "error: missing required -n / --input"
    usage
    exit 2
  fi
  if [[ ! -f "${nanopore}" ]] ; then
    >&2 echo -e "error: input '${nanopore}' does not exist"
    exit 2
  fi

  if [[ -z "${illumina}" ]]; then
    >&2 echo -e "error: missing required -i / --input"
    usage
    exit 2
  fi
  if [[ ! -f "${illumina}" ]] ; then
    >&2 echo -e "error: input '${illumina}' does not exist"
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
  local -r paramIllumina="$(realpath "${2}")"
  local -r paramOutput="$(realpath "${3}")"
  local -r paramProfile="${4}"

  $WORK_DIR/nextflow-23.10.0-all PhasingPipe.nf -c nextflow.config -profile $paramProfile --nanopore $paramNanopore --illumina $paramIllumina --output $paramOutput
}

main() {
  local nanopore=""
  local illumina=""
  local output=""
  local profile=""
  if command -v sbatch &> /dev/null; then
    profile="slurm"
  else
    profile="local"
  fi
  local args=$(getopt -o n:i:o: --long nanopore:,illumina:,output:,help: -- "$@")
  eval set -- "${args}"

  while :; do
    case "$1" in
    -n | --nanopore)
      nanopore="$2"
      shift 2 ;;
    -i | --illumina )
      illumina="$2"
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

  validate "${nanopore}" "${illumina}" "${output}"

  if ! [[ -d "${output}" ]]; then
    mkdir -p "${output}"
  fi

  execute_workflow "${nanopore}" "${illumina}" "${output}" "${profile}"
}

main "${@}"
