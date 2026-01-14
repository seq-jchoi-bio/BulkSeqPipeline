#!/bin/bash
# Converting BAM to BED (UCSC style) for ChIP-seq (local run)

set -o errexit
set -o nounset
set -o pipefail

# -----------------------
# Inputs
# -----------------------
BAMF="${1:?Error: missing input BAM (*.bam)}"
genome="${2:?Error: missing genome prefix (e.g., os*, hs*, ms*)}"
CPU_IN=${3:-2}

if [ ! -f "$BAMF" ]; then
    echo "Error: Input BAM does not exist: $BAMF"
    exit 1
fi

# -----------------------
# Variables
# -----------------------
bedgz="${BAMF%.bam}.bed.gz"
SAMCPUS="$CPU_IN"

# -----------------------
# Checks
# -----------------------
if ! command -v bamToBed >/dev/null 2>&1; then
    echo "Error: bamToBed not found in PATH (bedtools required)."
    exit 1
fi

echo "=========================================="
echo "Converting BAM -> BED.GZ"
echo "  BAM     : $BAMF"
echo "  Genome  : $genome"
echo "  BED.GZ  : $bedgz"
echo "  Threads : $SAMCPUS"
echo "=========================================="

# -----------------------
# Run
# -----------------------
if [[ "$genome" == os* ]]; then
    # NOTE: Use name-sorted stream for bamToBed
    samtools sort -@ "$SAMCPUS" -n -m 2G -O BAM "$BAMF" \
    | bamToBed -i - \
    | awk -v OFS="\t" '
    BEGIN {
        main_chromosomes["1"]; main_chromosomes["2"]; main_chromosomes["3"];
        main_chromosomes["4"]; main_chromosomes["5"]; main_chromosomes["6"];
        main_chromosomes["7"]; main_chromosomes["8"]; main_chromosomes["9"];
        main_chromosomes["10"]; main_chromosomes["11"]; main_chromosomes["12"];
        main_chromosomes["Mt"]; main_chromosomes["Pt"];
    }
    {
        if ($1 in main_chromosomes) {
            $1 = "Chr" $1;
        } else if ($1 ~ /^AP/ || $1 ~ /^AC/) {
            $1 = "ChrUn";
        } else if ($1 ~ /^Syng/) {
            $1 = "ChrSy";
        }

        if ($1 == "ChrMt" || $1 == "ChrPt" || $1 == "ChrUn" || $1 == "ChrSy") {
            next;
        }

        $4 = "read" NR;
        print $0;
    }' \
    | gzip -c > "$bedgz"

elif [[ "$genome" == hs* || "$genome" == ms* ]]; then
    samtools sort -@ "$SAMCPUS" -n -m 2G -O BAM "$BAMF" \
    | bamToBed -i - \
    | awk -v OFS="\t" '
        $1 !~ /^chrM$/ && $1 !~ /_/ { print $0; }
    ' \
    | gzip -c > "$bedgz"

else
    echo "Unsupported genome type: $genome"
    exit 1
fi

echo "Conversion complete: $bedgz"
