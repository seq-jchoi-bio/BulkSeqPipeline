#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

start_time=$(date "+%a %b %d %I:%M:%S %p %Z %Y")
VERSION=2026.01

HEADER="\n ==========================================="
HEADER+="\n   A Standard ChIP-Seq Pipeline"
HEADER+="\n   Version ${VERSION}"
HEADER+="\n   Pipeline created by: Janghyun Choi"
HEADER+="\n   Contact: jchoi@inha.ac.kr"
HEADER+="\n ==========================================="
echo -e "$HEADER"

printf "\n"
read -e -p "   1. Enter path to forward file (*.fastq.gz): " file1
printf "\n"
read -e -p "   2. Enter path to reverse file (*.fastq.gz): " file2
printf "\n"
read -e -p "   3. Enter genome prefix name (no path, auto-detection): " genome
printf "\n"
read -e -p "   4. Enter experiment name (no path): " expname
printf "\n"

base_dir="$(pwd)"

while true; do
    read -p "   5. Select assay type ('histone' or 'TF'): " assaytype
    case "$assaytype" in
        histone|TF ) break;;
        * ) echo "Invalid choice. Please enter 'histone' or 'TF'.";;
    esac
done

echo "==========================================="
echo "   System resource status"
echo "==========================================="

if command -v nproc >/dev/null 2>&1; then
  cpu_total="$(nproc)"
else
  cpu_total="$(getconf _NPROCESSORS_ONLN)"
fi
echo "   - CPU threads available: ${cpu_total}"

if command -v free >/dev/null 2>&1; then
  echo "   - Memory (free -h):"
  free -h
else
  if command -v sysctl >/dev/null 2>&1; then
    mem_bytes="$(sysctl -n hw.memsize 2>/dev/null || echo 0)"
    if [[ "$mem_bytes" -gt 0 ]]; then
      mem_gb="$(( mem_bytes / 1024 / 1024 / 1024 ))"
      echo "   - Memory total: ~${mem_gb} GB"
    fi
  fi
fi

echo "==========================================="
printf "\n"

while true; do
  read -p "   6. How many CPU cores/threads to use? (1-${cpu_total}): " cpu
  if [[ "$cpu" =~ ^[0-9]+$ ]] && [[ "$cpu" -ge 1 ]] && [[ "$cpu" -le "$cpu_total" ]]; then
    break
  fi
  echo "   Please enter a number between 1 and ${cpu_total}."
done
printf "\n"

# Resolve repository root (main script is placed at repo root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR}"

# Check input FASTQs
if [ ! -f "$file1" ]; then
    echo "Error: Forward file '$file1' does not exist."
    exit 1
fi
if [ ! -f "$file2" ]; then
    echo "Error: Reverse file '$file2' does not exist."
    exit 1
fi

# Check peak detection
while true; do
    echo -e "\n   ==========================================="
    echo "   [NOTE] Select 'no' if this is the input sample."
    echo "   [NOTE] Select 'no' if this is the target sample but do not want to perform peak-call."
    echo "   ==========================================="
    read -p "   Do you want to perform peak detection (MACS2)? (yes/no): " peak_choice
    case "$peak_choice" in
        [Yy]* ) peak_detection="yes"; break;;
        [Nn]* ) peak_detection="no"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

# Peak call logic (input control)
inputBed=""
inputExp=""
inputPath=""
if [ "$peak_detection" == "yes" ]; then
    echo -e "\n   ==========================================="
    echo "   [Note] Please enter the folder name (expname) containing the input file, including the full path."
    echo "   [Note] Ensure that the path includes only the folder name where the input file is located."
    echo "   [Note] This script recognizes only the file format: (expname)/(expname).flt.rd.bed.gz"
    echo "   ==========================================="

    while true; do
        read -e -p "   7. Enter the folder name (expname) containing the input file (relative path only): " inputBed

        if [ -z "$inputBed" ]; then
            echo "   Error: The folder (expname) cannot be empty. Please enter a valid name."
            continue
        fi

        inputBed="$(pwd)/$inputBed"
        inputExp="$(basename "$inputBed")"
        inputPath="$inputBed/$inputExp.flt.rd.bed.gz"

        if [ ! -d "$inputBed" ]; then
            echo "   Error: The directory '$inputBed' does not exist. Please enter a valid path."
            continue
        fi

        if [ -f "$inputPath" ]; then
            echo "   Found file: $inputPath"
            break
        else
            echo "   Error: The file '$inputPath' does not exist. Please enter a valid expname."
        fi
    done
fi

outdir="${base_dir}/${expname}"
mkdir -p "$outdir"

# Log directory
logdir="${outdir}/log"
mkdir -p "$logdir"

SCRIPTDIR="${REPO_ROOT}/bin"

bowtieidx="${REPO_ROOT}/refGenome/bowtie2_index/${genome}/genome"

trimdir="$outdir/Trim"
mkdir -p "$trimdir"
trimmed1="$trimdir/Trim_$(basename "${file1%.fastq.gz}").fastq.gz"
trimmed2="$trimdir/Trim_$(basename "${file2%.fastq.gz}").fastq.gz"
bamfile="${outdir}/${expname}.bam"
fbamfile="${outdir}/${expname}.flt.rd.bam"
bedfile="${outdir}/${expname}.flt.rd.bed.gz"

echo -e "\n **Output has been created in $outdir, and all results will be stored there.**"
echo -e "\n"
echo "==========================================="
echo "     - Trimming results, including QC, will be stored in: $trimdir"
echo -e "\n     - BAM files will be stored in: $outdir"
echo -e "\n     - BAM QC files will be stored in: ${outdir}_QC"
if [ "$peak_detection" == "yes" ]; then
    echo -e "\n     - Peak-call results will be stored in: ${outdir}/macs2_${expname}"
fi
echo -e "\n     - Pipeline log files will be stored in: ${outdir}/log"
echo -e "\n     - The final report will be generated in: ${outdir}/pipeline_summary.txt"
echo -e "\n     - Threads to use: ${cpu}"
echo "==========================================="

# 1. Trimming and QC (local)
echo -e "\n =========================================="
echo "   Running: trimmomaticChip.sh"
echo "=========================================="
bash "${SCRIPTDIR}/trimmomaticChip.sh" "${file1}" "${file2}" "${trimdir}" "${cpu}" \
  2>&1 | tee "${logdir}/01_trimmomatic.log"
echo -e "\n =========================================="
echo "   Done: trimmomaticChip.sh"
echo "=========================================="

# 2. Mapping (local)
echo -e "\n =========================================="
echo "   Running: bowtie2Chip.sh"
echo "=========================================="
bash "${SCRIPTDIR}/bowtie2Chip.sh" "${trimmed1}" "${trimmed2}" "${bowtieidx}" "${bamfile}" "${assaytype}" "${cpu}" \
  2>&1 | tee "${logdir}/02_bowtie2.log"
echo -e "\n =========================================="
echo "   Done: bowtie2Chip.sh"
echo "=========================================="

# 3. Converting/QC (local)
echo -e "\n =========================================="
echo "   Running: filterQCChip.sh"
echo "=========================================="
bash "${SCRIPTDIR}/filterQCChip.sh" "${bamfile}" "${cpu}" \
  2>&1 | tee "${logdir}/03_filterQC.log"
echo -e "\n =========================================="
echo "   Done: filterQCChip.sh"
echo "=========================================="

# 4. BED (local)
echo -e "\n =========================================="
echo "   Running: convertingBEDChip.sh"
echo "=========================================="
bash "${SCRIPTDIR}/convertingBEDChip.sh" "${fbamfile}" "${genome}" "${cpu}" \
  2>&1 | tee "${logdir}/04_bed.log"
echo -e "\n =========================================="
echo "   Done: convertingBEDChip.sh"
echo "=========================================="

# 5. MACS2 (conditional, local)
if [ "$peak_detection" == "yes" ]; then
    echo -e "\n =========================================="
    echo "   Running: macs2Chip.sh"
    echo "=========================================="
    bash "${SCRIPTDIR}/macs2Chip.sh" "${genome}" "${expname}" "${bedfile}" "${inputPath}" "${assaytype}" \
      2>&1 | tee "${logdir}/05_macs2.log"
    echo -e "\n =========================================="
    echo "   Done: macs2Chip.sh"
    echo "=========================================="
else
    echo -e "\n   Skipping MACS2 step (peak_detection=no)."
fi

# 6. Generate report (local)
echo -e "\n =========================================="
echo "   Running: generate_reportChip.sh"
echo "=========================================="
bash "${SCRIPTDIR}/generate_reportChip.sh" "${expname}" "${start_time}" \
  2>&1 | tee "${logdir}/06_report.log"
echo -e "\n =========================================="
echo "   Done: generate_reportChip.sh"
echo "=========================================="
