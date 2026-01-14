#!/bin/bash
# Mapping using HISAT2 for RNA-seq (local execution)

set -o errexit
set -o nounset

# Arguments
SEQF1=$1        # Trimmed R1 fastq.gz
SEQF2=$2        # Trimmed R2 fastq.gz
IDXPREFIX=$3    # HISAT2 index prefix
OUTBAMF=$4      # Output BAM file
CPU_IN=${5:-2}

# Threads (adjustable)
NCPUS="$CPU_IN"
SAMCPUS="$CPU_IN"

# Check inputs
if [ ! -f "$SEQF1" ]; then
    echo "Error: Input file $SEQF1 does not exist."
    exit 1
fi

if [ ! -f "$SEQF2" ]; then
    echo "Error: Input file $SEQF2 does not exist."
    exit 1
fi

if [ ! -f "${IDXPREFIX}.1.ht2" ]; then
    echo "Error: HISAT2 index prefix '$IDXPREFIX' not found."
    exit 1
fi

# Output directory
mkdir -p "$(dirname "$OUTBAMF")"

# Run HISAT2 + BAM conversion
hisat2 \
    -p $NCPUS \
    -x "$IDXPREFIX" \
    -1 "$SEQF1" \
    -2 "$SEQF2" \
| samtools view -@ $SAMCPUS -b - > "$OUTBAMF"

echo "HISAT2 mapping finished: $OUTBAMF"
