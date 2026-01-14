#!/usr/bin/env bash
set -euo pipefail

EFFECTIVE_GENOME_SIZE=373245519
DONE_WORD="done"

trap 'echo "[ERROR] Failed at line $LINENO: $BASH_COMMAND" >&2' ERR
die() { echo "[ERROR] $*" >&2; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }

# Global arrays (no eval, no nameref)
BAM_ITEMS=()
REGION_ITEMS=()

trim_cr() {
  # Remove possible CR from Windows-style input
  printf '%s' "$1" | sed 's/\r$//'
}

lower() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

default_out_bw() {
  local bam="$1"
  local base
  base="$(basename "$bam")"
  base="${base%.bam}"
  echo "${base}.bigwig"
}

default_out_bw_region() {
  local bam="$1"
  local region="$2"
  local base rtag
  base="$(basename "$bam")"
  base="${base%.bam}"
  rtag="$(echo "$region" | sed 's/:/_/g')"
  echo "${base}_${rtag}.bigwig"
}

resolve_bam_from_item() {
  # folder -> ./folder/folder.flt.rd.bam
  # file   -> use directly
  local item="$1"
  if [[ -d "$item" ]]; then
    local b
    b="$(basename "$item")"
    echo "${item}/${b}.flt.rd.bam"
  else
    echo "$item"
  fi
}

detect_chr_prefix() {
  local bam="$1"
  local sn
  sn="$(samtools view -H "$bam" | awk -F'\t' '$1=="@SQ"{for(i=1;i<=NF;i++) if($i~/^SN:/){sub("SN:","",$i); print $i; exit}}')"
  [[ -z "$sn" ]] && die "Failed to read @SQ SN from BAM header: $bam"
  [[ "$sn" == Chr* ]] && echo "Chr" || echo ""
}

read_bam_items_until_done() {
  BAM_ITEMS=()
  echo "   Type '${DONE_WORD}' to finish input."
  while true; do
    read -r line || break
    line="$(trim_cr "$line")"
    [[ -z "$line" ]] && continue
    [[ "$(lower "$line")" == "$DONE_WORD" ]] && break
    BAM_ITEMS+=("$line")
  done
  [[ "${#BAM_ITEMS[@]}" -ge 1 ]] || die "No BAM items provided before '${DONE_WORD}'."
}

read_regions_until_done() {
  REGION_ITEMS=()
  echo "   Type '${DONE_WORD}' to finish input."
  while true; do
    read -r line || break
    line="$(trim_cr "$line")"
    [[ -z "$line" ]] && continue
    [[ "$(lower "$line")" == "$DONE_WORD" ]] && break
    REGION_ITEMS+=("$line")
  done
  [[ "${#REGION_ITEMS[@]}" -ge 1 ]] || die "No regions provided before '${DONE_WORD}'."
}

need_cmd samtools
need_cmd bamCoverage

echo "==========================================="
echo " bigWig generator (bamCoverage)"
echo " - low:  bs=100"
echo " - high: bs=10 + smoothLength=20 (very slow)"
echo "==========================================="

# Resolution
while true; do
  read -r -p "1) Choose resolution [low/high]: " RES
  RES="$(lower "$(trim_cr "$RES")")"
  [[ "$RES" == "low" || "$RES" == "high" ]] && break
done

BS=100
SMOOTH_ARGS=()
if [[ "$RES" == "high" ]]; then
  BS=10
  SMOOTH_ARGS=(--smoothLength 20)
  echo "!! WARNING: High-resolution mode is very slow."
fi

# Mode
while true; do
  read -r -p "2) Track mode [all/specific]: " MODE
  MODE="$(lower "$(trim_cr "$MODE")")"
  [[ "$MODE" == "all" || "$MODE" == "specific" ]] && break
done

# Cores
while true; do
  read -r -p "3) Number of cores (required): " CORE
  CORE="$(trim_cr "$CORE")"
  [[ "$CORE" =~ ^[0-9]+$ ]] && [[ "$CORE" -ge 1 ]] && break
done

if [[ "$MODE" == "all" ]]; then
  echo "4) Enter sample folders or BAM files (one per line)"
  echo "   - folder  -> ./folder/folder.flt.rd.bam"
  echo "   - bamfile -> used directly"
  echo
  read_bam_items_until_done

  for ITEM in "${BAM_ITEMS[@]}"; do
    BAM="$(resolve_bam_from_item "$ITEM")"
    echo "[INFO] Resolved: $ITEM -> $BAM"
    [[ -f "$BAM" ]] || die "BAM not found: $BAM"

    OUT="$(default_out_bw "$BAM")"
    echo "[RUN] bamCoverage -b \"$BAM\" -o \"$OUT\" ..."
    bamCoverage -b "$BAM" -o "$OUT" \
      -bs "$BS" --normalizeUsing RPGC -e \
      --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
      -p "$CORE" "${SMOOTH_ARGS[@]}"
  done

else
  while true; do
    read -r -p "4) Enter sample folder or BAM file: " ITEM
    ITEM="$(trim_cr "$ITEM")"
    [[ -z "$ITEM" ]] && continue
    BAM="$(resolve_bam_from_item "$ITEM")"
    echo "[INFO] Resolved: $ITEM -> $BAM"
    [[ -f "$BAM" ]] && break
    echo "   BAM not found: $BAM"
  done

  PREFIX="$(detect_chr_prefix "$BAM")"
  if [[ -n "$PREFIX" ]]; then
    echo "5) BAM header uses 'Chr' prefix. Example: Chr1:1000:2000"
  else
    echo "5) BAM header has no 'Chr' prefix. Example: 1:1000:2000"
  fi

  echo "6) Enter regions (one per line)"
  echo
  read_regions_until_done

  for REGION in "${REGION_ITEMS[@]}"; do
    # Basic format validation
    if [[ -n "$PREFIX" ]]; then
      [[ "$REGION" == Chr*:*:* ]] || die "Region must start with 'Chr' (e.g., Chr1:100:200): $REGION"
    else
      [[ "$REGION" != Chr*:*:* ]] || die "Region must NOT start with 'Chr' (e.g., 1:100:200): $REGION"
      [[ "$REGION" == *:*:* ]] || die "Region must be like '1:100:200': $REGION"
    fi

    OUTDIR="$(echo "$REGION" | sed 's/:/_/g')"
    mkdir -p "$OUTDIR"
    OUT="${OUTDIR}/$(default_out_bw_region "$BAM" "$REGION")"

    echo "[RUN] bamCoverage -b \"$BAM\" -o \"$OUT\" --region \"$REGION\" ..."
    bamCoverage -b "$BAM" -o "$OUT" \
      --region "$REGION" \
      -bs "$BS" --normalizeUsing RPGC -e \
      --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" \
      -p "$CORE" "${SMOOTH_ARGS[@]}"
  done
fi

echo "Done."
