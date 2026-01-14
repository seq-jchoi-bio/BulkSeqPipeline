#!/bin/bash

set -o errexit
set -o nounset
set -o pipefail

start_time=$(date "+%a %b %d %I:%M:%S %p %Z %Y")
VERSION=2026.01

HEADER="\n ==========================================="
HEADER+="\n   A Standard Bulk RNA-Seq Pipeline"
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
read -e -p "   5. Enter GTF file (no path, *.gtf, auto-detection): " gtffile
printf "\n"

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

base_dir="$(pwd)"

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

# Check gtf file
gtfpath="${REPO_ROOT}/refGenome/gtf/${genome}/${gtffile}"
if [ ! -f "$gtfpath" ]; then
    echo "Error: GTF file '$gtfpath' does not exist."
    exit 1
fi

# Counting logic
while true; do
    printf "\n"
    read -p "   Do you want to perform counting (htseq)? (yes/no): " count_choice
    case "$count_choice" in
        [Yy]* ) counting="yes"; break;;
        [Nn]* ) counting="no"; break;;
        * ) echo "Please answer yes or no.";;
    esac
done

# Counting condition
htseqBam=()
if [ "$counting" == "yes" ]; then
    printf "\n"
    echo "   Enter *BAM file expriment names* (press Enter after each, type 'done' to finish):"
    while true; do
        read -e -p "   > " bamfile
        if [ "$bamfile" == "done" ]; then
            break
        fi
        htseqBam+=("$bamfile")
    done

    mkdir -p "${base_dir}/Count"
fi

outdir="${base_dir}/${expname}"
mkdir -p "$outdir"

# Log directory
logdir="${outdir}/log"
mkdir -p "$logdir"

if [ "$counting" == "yes" ]; then
    echo "Count result directory created: ${base_dir}/Count"
fi

SCRIPTDIR="${REPO_ROOT}/bin"
hisat2idx="${REPO_ROOT}/refGenome/hisat2_index/${genome}/genome"

trimdir="$outdir/Trim"
mkdir -p "$trimdir"
trimmed1="$trimdir/Trim_$(basename ${file1%.fastq.gz}).fastq.gz"
trimmed2="$trimdir/Trim_$(basename ${file2%.fastq.gz}).fastq.gz"
bamfile="${outdir}/${expname}.bam"

echo -e "\n **Output has been created in $outdir, and all results will be stored there.**"
echo -e "\n"
echo "==========================================="
echo -e "   - Trimming results, including QC, will be stored in: $trimdir"
echo -e "\n   - BAM files will be stored in: $outdir"
echo -e "\n   - BAM QC files will be stored in: ${outdir}/${expname}_QC"
if [ "$counting" == "yes" ]; then
    echo -e "\n   - Count result will be generated in: ${base_dir}/Count/count.txt"
fi
echo -e "\n   - Pipeline log files will be stored in: ${outdir}/log"
echo -e "\n   - The final report will be generated in: ${outdir}/pipeline_summary.txt"
echo -e "\n   - Threads to use: ${cpu}"
echo "=========================================="

# Trimming and QC (local)
echo -e "\n =========================================="
echo "   Running: trimmomatic_rna.sh"
echo "=========================================="
bash "${SCRIPTDIR}/trimmomatic_rna.sh" "${file1}" "${file2}" "${trimdir}" "${cpu}" \
  2>&1 | tee "${logdir}/01_trimmomatic.log"
echo -e "\n =========================================="
echo "   Done: trimmomatic_rna.sh"
echo "=========================================="

# Mapping (local)
echo -e "\n =========================================="
echo "   Running: hisat2_rna.sh"
echo "=========================================="
bash "${SCRIPTDIR}/hisat2_rna.sh" "${trimmed1}" "${trimmed2}" "${hisat2idx}" "${bamfile}" "${cpu}" \
  2>&1 | tee "${logdir}/02_hisat2.log"
echo -e "\n =========================================="
echo "   Done: hisat2_rna.sh"
echo "=========================================="

# Converting/QC (local)
echo -e "\n =========================================="
echo "   Running: filterQC_rna.sh"
echo "=========================================="
bash "${SCRIPTDIR}/filterQC_rna.sh" "${bamfile}" "${cpu}" \
  2>&1 | tee "${logdir}/03_filterQC.log"
echo -e "\n =========================================="
echo "   Done: filterQC_rna.sh"
echo "=========================================="

# Counting (optional, local)
if [[ "$counting" == "yes" && ${#htseqBam[@]} -gt 0 ]]; then
    echo "=========================================="
    echo "   Running htseq_rna.sh with inputs:"
    for value in "${htseqBam[@]}"; do
        echo "   ${base_dir}/${value}/${value}.flt.rd.bam"
    done
    echo "=========================================="

    echo -e "\n =========================================="
    echo "   Running: htseq_rna.sh"
    echo "=========================================="
    bash "${SCRIPTDIR}/htseq_rna.sh" "${base_dir}/Count" "$gtfpath" "${htseqBam[@]}" \
      2>&1 | tee "${logdir}/04_htseq.log"
    echo -e "\n =========================================="
    echo "   Done: htseq_rna.sh"
    echo "=========================================="
fi

# Generate report (local)
echo -e "\n =========================================="
echo "   Running: generate_report.sh"
echo "=========================================="
bash "${SCRIPTDIR}/generate_report.sh" "${expname}" "${start_time}" \
  2>&1 | tee "${logdir}/05_report.log"
echo -e "\n =========================================="
echo "   Done: generate_report.sh"
echo "=========================================="
