#!/bin/bash
# Peak calling via MACS2 for ChIP-seq (local run)

set -o errexit
set -o nounset
set -o pipefail

# -----------------------
# Inputs
# -----------------------
GENOME="${1:?Error: missing genome prefix (e.g., os*, hs*, mm*)}"
OUTPREFIX="${2:?Error: missing experiment name (outprefix)}"
TargetBED="${3:?Error: missing target BED (.bed.gz)}"
INPUTBED="${4:?Error: missing input BED (.bed.gz)}"
AssayType="${5:?Error: missing assay type (TF|histone)}"

# -----------------------
# Output dir (match main pipeline expectation)
# main pipeline prints: ${outdir}/macs2_${expname}
# and calls: macs2Chip.sh ${genome} ${expname} ${bedfile} ${inputPath} ${assaytype}
# so we place results under: $PWD/${OUTPREFIX}/macs2_${OUTPREFIX}
# -----------------------
OUTDIR="$PWD/${OUTPREFIX}"
MACS2DIR="${OUTDIR}/macs2_${OUTPREFIX}"
mkdir -p "$MACS2DIR"

# -----------------------
# Checks
# -----------------------
if [ ! -f "$TargetBED" ]; then
    echo "Error: Target BED not found: $TargetBED"
    exit 1
fi
if [ ! -f "$INPUTBED" ]; then
    echo "Error: Input BED not found: $INPUTBED"
    exit 1
fi

# -----------------------
# Genome size mapping
# -----------------------
gsize="${GENOME:0:2}"
if [[ "$gsize" == "hg" ]]; then
    gsize="hs"
elif [[ "$gsize" == "mm" ]]; then
    gsize="mm"
elif [[ "$gsize" == "os" ]]; then
    gsize="373245519"
fi

echo "=========================================="
echo "Running MACS2 callpeak"
echo "  Genome    : $GENOME  (gsize=$gsize)"
echo "  AssayType : $AssayType"
echo "  Name      : $OUTPREFIX"
echo "  TargetBED : $TargetBED"
echo "  InputBED  : $INPUTBED"
echo "  Outdir    : $MACS2DIR"
echo "=========================================="

# -----------------------
# Run
# -----------------------
if [[ "$AssayType" == "TF" ]]; then
    macs2 callpeak -n "$OUTPREFIX" \
        -g "$gsize" --nomodel \
        --shift -75 --extsize 150 \
        -f BED -B --SPMR \
        -t "$TargetBED" -c "$INPUTBED" \
        --outdir "$MACS2DIR"

elif [[ "$AssayType" == "histone" ]]; then
    macs2 callpeak -n "$OUTPREFIX" \
        -g "$gsize" --nomodel \
        --shift 0 --extsize 147 \
        --broad --broad-cutoff 0.1 \
        -f BED -B --SPMR \
        -t "$TargetBED" -c "$INPUTBED" \
        --outdir "$MACS2DIR"

else
    echo "Invalid AssayType: $AssayType. Use 'TF' or 'histone'."
    exit 1
fi

# -----------------------
# Compress bdg outputs (if present)
# -----------------------
TREAT_BDG="${MACS2DIR}/${OUTPREFIX}_treat_pileup.bdg"
CTRL_BDG="${MACS2DIR}/${OUTPREFIX}_control_lambda.bdg"

if [ -f "$TREAT_BDG" ]; then
    gzip -f "$TREAT_BDG"
fi
if [ -f "$CTRL_BDG" ]; then
    gzip -f "$CTRL_BDG"
fi

echo "Done."
