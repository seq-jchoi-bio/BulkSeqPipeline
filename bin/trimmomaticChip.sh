#!/bin/bash
# Trimming using Trimmomatic for ChIP-seq (local run)

set -o errexit
set -o nounset
set -o pipefail

# -----------------------
# Inputs
# -----------------------
PAIRED_IN1="${1:?Error: missing R1 fastq.gz}"
PAIRED_IN2="${2:?Error: missing R2 fastq.gz}"
OUTDIR="${3:?Error: missing output directory}"
CPU_IN=${4:-2}

# -----------------------
# Resolve repo root (bin/ script assumed)
# -----------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TRIMMOMATIC_JAR="${REPO_ROOT}/programs/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERLC="${REPO_ROOT}/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# -----------------------
# Parameters
# -----------------------
CUT="2:30:10"
LEAD="3"
TRAIL="3"
WINDOW="4:15"
LENGTH="36"
TCORE="$CPU_IN"

# -----------------------
# Checks
# -----------------------
if [ ! -f "$PAIRED_IN1" ]; then
    echo "Error: Input file $PAIRED_IN1 does not exist."
    exit 1
fi
if [ ! -f "$PAIRED_IN2" ]; then
    echo "Error: Input file $PAIRED_IN2 does not exist."
    exit 1
fi

mkdir -p "$OUTDIR"

if [ ! -f "$TRIMMOMATIC_JAR" ]; then
    echo "Error: Trimmomatic jar not found: $TRIMMOMATIC_JAR"
    exit 1
fi
if [ ! -f "$ADAPTERLC" ]; then
    echo "Error: Adapter file not found: $ADAPTERLC"
    exit 1
fi

# -----------------------
# Outputs
# -----------------------
TRIM_OUT1="$OUTDIR/Trim_$(basename "${PAIRED_IN1%.fastq.gz}").fastq.gz"
TRIM_OUT2="$OUTDIR/Trim_$(basename "${PAIRED_IN2%.fastq.gz}").fastq.gz"
UNTRIM_OUT1="$OUTDIR/unTrim_$(basename "${PAIRED_IN1%.fastq.gz}").fastq.gz"
UNTRIM_OUT2="$OUTDIR/unTrim_$(basename "${PAIRED_IN2%.fastq.gz}").fastq.gz"

# -----------------------
# Run
# -----------------------
echo "=========================================="
echo "Running Trimmomatic (PE)"
echo "  REPO_ROOT : $REPO_ROOT"
echo "  R1        : $PAIRED_IN1"
echo "  R2        : $PAIRED_IN2"
echo "  OUTDIR    : $OUTDIR"
echo "  Threads   : $TCORE"
echo "=========================================="

java -jar "$TRIMMOMATIC_JAR" \
    PE -threads "$TCORE" -phred33 \
    "$PAIRED_IN1" "$PAIRED_IN2" \
    "$TRIM_OUT1" "$UNTRIM_OUT1" "$TRIM_OUT2" "$UNTRIM_OUT2" \
    "ILLUMINACLIP:$ADAPTERLC:$CUT" \
    "LEADING:$LEAD" \
    "TRAILING:$TRAIL" \
    "SLIDINGWINDOW:$WINDOW" \
    "MINLEN:$LENGTH"

echo "=========================================="
echo "Running FastQC on trimmed reads"
echo "=========================================="

mkdir -p "$OUTDIR/trimmed_fastQC"
fastqc -o "$OUTDIR/trimmed_fastQC" -t "$TCORE" "$TRIM_OUT1" "$TRIM_OUT2"

echo "Done."
