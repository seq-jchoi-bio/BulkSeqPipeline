#!/bin/bash
# Mapping to reference genome using Bowtie2 for ChIP-seq (local run)

set -o errexit
set -o nounset
set -o pipefail

# -----------------------
# Inputs
# -----------------------
SEQF1="${1:?Error: missing trimmed R1 fastq.gz}"
SEQF2="${2:?Error: missing trimmed R2 fastq.gz}"
IDXPREFIX="${3:?Error: missing bowtie2 index prefix}"
OUTBAMF="${4:?Error: missing output BAM file}"
ArrayType="${5:?Error: missing assay type (histone|TF)}"
CPU_IN=${6:-2}

# -----------------------
# Resources
# -----------------------
NCPUS="$CPU_IN"
SAMCPUS="$CPU_IN"

# -----------------------
# Checks
# -----------------------
if [ ! -f "$SEQF1" ]; then
    echo "Error: File not found: $SEQF1"
    exit 1
fi
if [ ! -f "$SEQF2" ]; then
    echo "Error: File not found: $SEQF2"
    exit 1
fi
if [ ! -f "${IDXPREFIX}.1.bt2" ] && [ ! -f "${IDXPREFIX}.1.bt2l" ]; then
    echo "Error: Bowtie2 index not found: ${IDXPREFIX}.*.bt2"
    exit 1
fi

OUTDIR="$(dirname "$OUTBAMF")"
mkdir -p "$OUTDIR"

# -----------------------
# Run
# -----------------------
echo "=========================================="
echo "Running Bowtie2 (ChIP-seq)"
echo "  Assay type : $ArrayType"
echo "  R1         : $SEQF1"
echo "  R2         : $SEQF2"
echo "  Index      : $IDXPREFIX"
echo "  Output BAM : $OUTBAMF"
echo "  Threads    : $NCPUS"
echo "=========================================="

case "$ArrayType" in
    histone)
        echo "Mode: --sensitive (histone)"
        bowtie2 -p "$NCPUS" --sensitive -x "$IDXPREFIX" \
            -1 "$SEQF1" -2 "$SEQF2" \
        | samtools view -@ "$SAMCPUS" -b -o "$OUTBAMF"
        ;;
    TF)
        echo "Mode: --very-sensitive (TF)"
        bowtie2 -p "$NCPUS" --very-sensitive -x "$IDXPREFIX" \
            -1 "$SEQF1" -2 "$SEQF2" \
        | samtools view -@ "$SAMCPUS" -b -o "$OUTBAMF"
        ;;
    *)
        echo "Error: Unknown assay type '$ArrayType' (use 'histone' or 'TF')."
        exit 1
        ;;
esac

echo "Done."
