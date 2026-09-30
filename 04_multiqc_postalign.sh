#!/usr/bin/env bash
set -euo pipefail

command -v multiqc >/dev/null 2>&1 || {
  echo "ERROR: multiqc not found on PATH" >&2
  exit 127
}

PASSED_LIST="qc_summary/passed_samples.txt"
FAILED_LIST="qc_summary/failed_samples.txt"
OUT_DIR="qc_reports_postalign"
LINKS_DIR="$OUT_DIR/_links"

[[ -s "$PASSED_LIST" || -s "$FAILED_LIST" ]] || {
  echo "ERROR: run QC classification first" >&2
  exit 2
}

mkdir -p "$LINKS_DIR/all" "$LINKS_DIR/passed" "$LINKS_DIR/failed"

link_qc_for_bam() {
  local bam="$1"
  local dest="$2"
  local base
  base="$(basename "$bam" .bam)"

  [[ -s "${base}_flagstat.txt" ]] &&
    ln -sf "$(readlink -f "${base}_flagstat.txt")" "$dest/${base}_flagstat.txt"

  [[ -d "${base}_qualimap_report" ]] &&
    ln -sfn "$(readlink -f "${base}_qualimap_report")" "$dest/${base}_qualimap_report"
}

populate() {
  local list="$1"
  local dest="$2"
  [[ -s "$list" ]] || return 0

  while IFS= read -r bam; do
    [[ -n "$bam" ]] || continue
    link_qc_for_bam "$bam" "$dest"
    link_qc_for_bam "$bam" "$LINKS_DIR/all"
  done < "$list"
}

populate "$PASSED_LIST" "$LINKS_DIR/passed"
populate "$FAILED_LIST" "$LINKS_DIR/failed"

[[ -n "$(find "$LINKS_DIR/all" -mindepth 1 -print -quit)" ]] &&
  multiqc "$LINKS_DIR/all" -o "$OUT_DIR" --filename multiqc_postalign_all

[[ -n "$(find "$LINKS_DIR/passed" -mindepth 1 -print -quit)" ]] &&
  multiqc "$LINKS_DIR/passed" -o "$OUT_DIR" --filename multiqc_postalign_passed

[[ -n "$(find "$LINKS_DIR/failed" -mindepth 1 -print -quit)" ]] &&
  multiqc "$LINKS_DIR/failed" -o "$OUT_DIR" --filename multiqc_postalign_failed
