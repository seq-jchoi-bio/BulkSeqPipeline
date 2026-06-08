#!/bin/bash
# Converting files and QC (local execution)

set -o errexit
set -o nounset

# Resolve repository root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
source "${REPO_ROOT}/bin/resource_utils.sh"

# Tools
QUALIMAP_BIN="${REPO_ROOT}/programs/qualimap_v2.3/qualimap"
PICARD_JAR="${REPO_ROOT}/programs/picard/picard.jar"

# Argument
BAMF=$1
CPU_IN=${2:-2}
SORT_MEM=${3:-}

# Parameters
prepare_resource_settings "$CPU_IN" "RNA filtering/QC" "$SORT_MEM"
SAMCPUS="$RESOURCE_SORT_THREADS"
SORT_MEM="$RESOURCE_SORT_MEM"
MAPQ=20

# Checks
if [ ! -f "$BAMF" ]; then
    echo "Error: Input BAM file $BAMF does not exist."
    exit 1
fi
if [ ! -f "$PICARD_JAR" ]; then
    echo "Error: Picard jar not found: $PICARD_JAR"
    exit 1
fi
if [ ! -x "$QUALIMAP_BIN" ]; then
    echo "Error: Qualimap binary not executable: $QUALIMAP_BIN"
    exit 1
fi

# Prefix
EXPNAME="${BAMF%.bam}"

# Output files
TBAMF="${EXPNAME}.tmp.bam"
IBAMF="${EXPNAME}.im.bam"
SBAMF="${EXPNAME}.flt.bam"
DBAMF="${EXPNAME}.flt.rd.bam"
METRICF="${EXPNAME}.flt.rd.metric.txt"
OUTQC="${EXPNAME}_QC"

# Filtering by MAPQ and coordinate sorting
echo "=========================================="
echo "Filtering (MAPQ >= ${MAPQ}) and sorting"
echo "  IN                : $BAMF"
echo "  Threads           : $SAMCPUS"
echo "  Sort memory/thread: $SORT_MEM"
echo "=========================================="

samtools view -@ "$SAMCPUS" -b -q "$MAPQ" "$BAMF" > "$TBAMF"
samtools sort -@ "$SAMCPUS" -m "$SORT_MEM" -O BAM -o "$IBAMF" "$TBAMF"

# Remove duplicates (Picard via jar)
java -jar "$PICARD_JAR" MarkDuplicates \
    I="$IBAMF" \
    O="$SBAMF" \
    M="$METRICF" \
    REMOVE_DUPLICATES=true

# Name-sort for counting
samtools sort -@ "$SAMCPUS" -m "$SORT_MEM" -n -O BAM -o "$DBAMF" "$SBAMF"

# Collect BAM stats
samtools stats "$BAMF"  > "${EXPNAME}.samstats.txt"
samtools stats "$SBAMF" > "${SBAMF%.bam}.samstats.txt"
samtools stats "$DBAMF" > "${DBAMF%.bam}.samstats.txt"

# Qualimap QC (coordinate-sorted BAM only)
"$QUALIMAP_BIN" bamqc -bam "$SBAMF" -outdir "$OUTQC"

# Cleanup
rm -f "$TBAMF" "$IBAMF" "$METRICF"
find . -name "samtools.*.tmp.*.bam" -type f -exec rm -f {} +
