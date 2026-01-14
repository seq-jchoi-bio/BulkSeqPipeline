#!/bin/bash
# Converting files and QC for ChIP-seq (local run)

set -o errexit
set -o nounset
set -o pipefail

# -----------------------
# Resolve repository root
# -----------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Tools (repo-relative)
QUALIMAP_BIN="${REPO_ROOT}/programs/qualimap_v2.3/qualimap"
PICARD_JAR="${REPO_ROOT}/programs/picard/picard.jar"

# -----------------------
# Input
# -----------------------
BAMF="${1:?Error: missing input BAM (*.bam)}"

if [ ! -f "$BAMF" ]; then
    echo "Error: Input BAM does not exist: $BAMF"
    exit 1
fi

CPU_IN=${2:-2}

# -----------------------
# Variables
# -----------------------
EXPNAME="${BAMF%.bam}"
TBAMF="${EXPNAME}.tmp.bam"
SBAMF="${EXPNAME}.flt.bam"
DBAMF="${EXPNAME}.flt.rd.bam"
METRICF="${EXPNAME}.flt.rd.metric.txt"
OUTQC="${EXPNAME}_QC"

SAMCPUS="$CPU_IN"
MAPQ="${MAPQ_CUTOFF:-20}"

# -----------------------
# Tool checks
# -----------------------
if [ ! -x "$QUALIMAP_BIN" ]; then
    echo "Error: Qualimap not found or not executable: $QUALIMAP_BIN"
    exit 1
fi
if [ ! -f "$PICARD_JAR" ]; then
    echo "Error: Picard jar not found: $PICARD_JAR"
    exit 1
fi

# -----------------------
# Filtering & sorting
# -----------------------
echo "=========================================="
echo "Filtering (MAPQ >= ${MAPQ}) and sorting"
echo "  IN      : $BAMF"
echo "  TMP     : $TBAMF"
echo "  SORTED  : $SBAMF"
echo "  Threads : $SAMCPUS"
echo "=========================================="

samtools view -@ "$SAMCPUS" -b -q "$MAPQ" "$BAMF" > "$TBAMF"
samtools sort -@ "$SAMCPUS" -O BAM -o "$SBAMF" "$TBAMF"
samtools index "$SBAMF"

# -----------------------
# Remove duplicates (Picard)
# -----------------------
echo "=========================================="
echo "Removing duplicates (Picard MarkDuplicates)"
echo "  IN  : $SBAMF"
echo "  OUT : $DBAMF"
echo "=========================================="

java -jar "$PICARD_JAR" MarkDuplicates \
    I="$SBAMF" \
    O="$DBAMF" \
    M="$METRICF" \
    REMOVE_DUPLICATES=true

samtools index "$DBAMF"

# -----------------------
# Collect BAM stats
# -----------------------
echo "=========================================="
echo "Collecting samtools stats"
echo "=========================================="

samtools stats "$BAMF"  > "${EXPNAME}.samstats.txt"
samtools stats "$SBAMF" > "${SBAMF%.bam}.samstats.txt"
samtools stats "$DBAMF" > "${DBAMF%.bam}.samstats.txt"

# Cleanup
rm -f "$TBAMF" "$METRICF"

# -----------------------
# Qualimap QC
# -----------------------
echo "=========================================="
echo "Running Qualimap bamqc"
echo "  BAM : $DBAMF"
echo "  OUT : $OUTQC"
echo "=========================================="

mkdir -p "$OUTQC"
"$QUALIMAP_BIN" bamqc -bam "$DBAMF" -outdir "$OUTQC"

echo "Done."
