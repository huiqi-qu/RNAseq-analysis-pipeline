#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --bam-list BAM_LIST --gtf ANNOTATION_GTF --qualimap QUALIMAP_EXECUTABLE [--threads N]" >&2
  exit 1
}

THREADS=4
BAM_LIST=""
ANNOTATION_GTF=""
QUALIMAP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-list) BAM_LIST="$2"; shift 2 ;;
    --gtf) ANNOTATION_GTF="$2"; shift 2 ;;
    --qualimap) QUALIMAP="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) usage ;;
  esac
done

[[ -n "$BAM_LIST" && -s "$BAM_LIST" ]] || usage
[[ -n "$ANNOTATION_GTF" && -s "$ANNOTATION_GTF" ]] || usage
[[ -n "$QUALIMAP" && -x "$QUALIMAP" ]] || usage

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found" >&2; exit 127; }
command -v featureCounts >/dev/null 2>&1 || { echo "ERROR: featureCounts not found" >&2; exit 127; }

while IFS= read -r BAM; do
  [[ -n "$BAM" ]] || continue
  [[ -s "$BAM" ]] || { echo "WARNING: missing BAM: $BAM" >&2; continue; }

  BASE="$(basename "$BAM" _Aligned.sortedByCoord.out.bam)"

  if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
    samtools index -@ "$THREADS" "$BAM"
  fi

  samtools flagstat -@ "$THREADS" "$BAM" > "${BASE}_flagstat.txt"

  "$QUALIMAP" bamqc \
    -bam "$BAM" \
    -outdir "${BASE}_qualimap_report" \
    -nt "$THREADS" \
    --java-mem-size=8G

  featureCounts \
    -T "$THREADS" \
    -a "$ANNOTATION_GTF" \
    -t exon \
    -g gene_id \
    -s 0 \
    -p -B -C \
    -o "${BASE}_gene_counts.txt" \
    "$BAM"
done < "$BAM_LIST"
