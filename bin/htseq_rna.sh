#!/bin/bash
# htseq-count (local execution)

set -o errexit
set -o nounset
set -o pipefail

output_dir=$1
gtfpath=$2
shift 2
htseqBam=("$@")

# Check output dir / gtf
mkdir -p "$output_dir"
if [ ! -f "$gtfpath" ]; then
    echo "Error: GTF file does not exist: $gtfpath"
    exit 1
fi

# base_dir: parent of Count dir (matches main pipeline design)
base_dir="$(cd "$(dirname "$output_dir")" && pwd)"

final_output="${output_dir}/count.txt"
STRANDOPT="no"
MODEOPT="intersection-nonempty"
ORDEROPT="name"

inputfiles=()

echo "Processing input files (in order):"
for value in "${htseqBam[@]}"; do
    inputfile="${base_dir}/${value}/${value}.flt.rd.bam"
    if [ -f "$inputfile" ]; then
        echo "Adding file: $inputfile"
        inputfiles+=("$inputfile")
    else
        echo "Error: File does not exist: $inputfile"
        exit 1
    fi
done

echo "Final input files in order: ${inputfiles[*]}"
echo "GTF file used: $gtfpath"
echo "Output will be saved to: $final_output"

# Run htseq-count
htseq-count -f bam --strand="$STRANDOPT" --mode="$MODEOPT" -r "$ORDEROPT" \
    "${inputfiles[@]}" "$gtfpath" > "$final_output"

echo "Final output saved to $final_output"
