#!/usr/bin/env bash
set -euo pipefail

BAM_LIST="${1:-}"
[[ -n "$BAM_LIST" && -s "$BAM_LIST" ]] || {
  echo "Usage: $0 BAM_LIST" >&2
  exit 1
}

MIN_MAPPED_PERCENTAGE=80
MIN_MEAN_MAPPING_QUALITY=30

mkdir -p qc_summary
PASSED="qc_summary/passed_samples.txt"
FAILED="qc_summary/failed_samples.txt"
METRICS="qc_summary/qc_metrics.tsv"

: > "$PASSED"
: > "$FAILED"
printf "bam\tmapped_percentage\tmean_mapping_quality\tstatus\n" > "$METRICS"

while IFS= read -r BAM; do
  [[ -n "$BAM" ]] || continue

  BASE="$(basename "$BAM" .bam)"
  FLAGSTAT="${BASE}_flagstat.txt"
  QUALIMAP_REPORT="${BASE}_qualimap_report/genome_results.txt"

  [[ -s "$FLAGSTAT" ]] || { echo "ERROR: missing $FLAGSTAT" >&2; exit 2; }
  [[ -s "$QUALIMAP_REPORT" ]] || { echo "ERROR: missing $QUALIMAP_REPORT" >&2; exit 2; }

  MAPPED=$(grep -m 1 "mapped (" "$FLAGSTAT" | awk '{print $5}' | tr -d '()%')
  MAPQ=$(grep "mean mapping quality" "$QUALIMAP_REPORT" | awk '{print $4}')

  STATUS="PASS"
  awk -v x="$MAPPED" -v min="$MIN_MAPPED_PERCENTAGE" 'BEGIN {exit !(x < min)}' && STATUS="FAIL"
  awk -v x="$MAPQ" -v min="$MIN_MEAN_MAPPING_QUALITY" 'BEGIN {exit !(x < min)}' && STATUS="FAIL"

  printf "%s\t%s\t%s\t%s\n" "$BAM" "$MAPPED" "$MAPQ" "$STATUS" >> "$METRICS"

  if [[ "$STATUS" == "PASS" ]]; then
    echo "$BAM" >> "$PASSED"
  else
    echo "$BAM" >> "$FAILED"
  fi
done < "$BAM_LIST"
