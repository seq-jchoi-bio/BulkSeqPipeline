#!/bin/bash
# Generating summary report for ChIP-seq (local execution)

set -o errexit
set -o nounset
set -o pipefail

EXPNAME="${1:?Error: missing experiment name}"
start_T="${2:?Error: missing start time string}"

end_time=$(date "+%a %b %d %I:%M:%S %p %Z %Y")

BASE_DIR="$(pwd)"
OUTDIR="${BASE_DIR}/${EXPNAME}"
LOGDIR="${OUTDIR}/log"
LOGFILE="${OUTDIR}/pipeline_summary.txt"

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

LOGDATA="Pipeline Summary Report\n"
LOGDATA+="========================\n"
LOGDATA+="Pipeline: Standard ChIP-seq Pipeline (local)\n"
LOGDATA+="User: $(whoami)\n"
LOGDATA+="Start Time: $start_T\n"
LOGDATA+="End Time: $end_time\n"
LOGDATA+="Working Directory: $BASE_DIR\n"
LOGDATA+="Output Directory: $OUTDIR\n"
LOGDATA+="========================\n\n"

# (A) Preferred: collect logs from ${OUTDIR}/log/*.log if present
if ls "${LOGDIR}"/*.log >/dev/null 2>&1; then
    for f in $(ls -1 "${LOGDIR}"/*.log | sort); do
        step=$(basename "$f")
        LOGDATA+="# ${step}\n"
        LOGDATA+="----------------------------------------\n"
        LOGDATA+="$(cat "$f")\n"
        LOGDATA+="========================================\n\n"
    done
else
    LOGDATA+="No step logs found in: ${LOGDIR}\n"
    LOGDATA+="(Expected: *.log created by the main pipeline via tee)\n"
    LOGDATA+="========================================\n\n"
fi

# (B) Backward-compatible: collect legacy *.out/*.err in current dir if present
has_legacy=false
for file in "${BASE_DIR}"/*.err "${BASE_DIR}"/*.out; do
    if [[ -f "$file" ]]; then
        has_legacy=true
        break
    fi
done

if [[ "$has_legacy" == true ]]; then
    LOGDATA+="\nLegacy SLURM-style logs found in ${BASE_DIR} (*.out/*.err)\n"
    LOGDATA+="========================================\n\n"
    for f in $(ls -1 "${BASE_DIR}"/*.err "${BASE_DIR}"/*.out 2>/dev/null | sort); do
        LOGDATA+="# $(basename "$f")\n"
        LOGDATA+="----------------------------------------\n"
        LOGDATA+="$(cat "$f")\n"
        LOGDATA+="========================================\n\n"
    done

    mkdir -p "${LOGDIR}/legacy_slurm_logs"
    mv "${BASE_DIR}"/*.err "${BASE_DIR}"/*.out "${LOGDIR}/legacy_slurm_logs/" 2>/dev/null || true
fi

echo -e "$LOGDATA" > "$LOGFILE"
echo "The report has been completed: $LOGFILE"

# Optional cleanup (safe)
find "$BASE_DIR" -name "samtools.*.tmp.*.bam" -type f -exec rm -f {} + 2>/dev/null || true
