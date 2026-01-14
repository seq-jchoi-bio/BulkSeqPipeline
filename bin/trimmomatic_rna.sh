#!/bin/bash
#Trimming using trimmomatic for RNA-seq

# environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TRIMMOMATIC_JAR="${REPO_ROOT}/programs/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Variables
PAIRED_IN1=$1
PAIRED_IN2=$2
OUTDIR=$3
CPU_IN=${4:-2}
TCORE="$CPU_IN"

# Check IPTs exist
if [ ! -f "$PAIRED_IN1" ]; then
    echo "Error: Input file $PAIRED_IN1 does not exist."
    exit 1
fi

if [ ! -f "$PAIRED_IN2" ]; then
    echo "Error: Input file $PAIRED_IN2 does not exist."
    exit 1
fi

# Check ./Trim directory
if [ ! -d "$OUTDIR" ]; then
    echo "The directory does not exist. Creating folder: $OUTDIR."
    mkdir -p $OUTDIR
else
    echo "The folder already exists. Files will be saved in: $OUTDIR."
fi

TRIM_OUT1="$OUTDIR/Trim_$(basename ${PAIRED_IN1%.fastq.gz}).fastq.gz"
TRIM_OUT2="$OUTDIR/Trim_$(basename ${PAIRED_IN2%.fastq.gz}).fastq.gz"
UNTRIM_OUT1="$OUTDIR/unTrim_$(basename ${PAIRED_IN1%.fastq.gz}).fastq.gz"
UNTRIM_OUT2="$OUTDIR/unTrim_$(basename ${PAIRED_IN2%.fastq.gz}).fastq.gz"
ADAPTERLC="${REPO_ROOT}/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
CUT=2:30:10
LEAD=3
TRAIL=3
WINDOW=4:15
LENGTH=36
TCORE=10

#Running
java -jar "${TRIMMOMATIC_JAR}" \
    PE -threads $TCORE -phred33 "$PAIRED_IN1" "$PAIRED_IN2" \
    "${TRIM_OUT1}" "${UNTRIM_OUT1}" "${TRIM_OUT2}" "${UNTRIM_OUT2}" \
    ILLUMINACLIP:"$ADAPTERLC":$CUT LEADING:$LEAD TRAILING:$TRAIL \
    SLIDINGWINDOW:$WINDOW MINLEN:$LENGTH && \
mkdir -p $OUTDIR/trimmed_fastQC
fastqc -o "$OUTDIR/trimmed_fastQC" -t $TCORE ${TRIM_OUT1} ${TRIM_OUT2}

