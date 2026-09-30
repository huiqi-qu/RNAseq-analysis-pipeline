#!/usr/bin/env bash
set -euo pipefail

COUNT_LIST="${1:-}"
OUTPUT="${2:-combined_counts.csv}"

[[ -n "$COUNT_LIST" && -s "$COUNT_LIST" ]] || {
  echo "Usage: $0 GENE_COUNT_FILE_LIST [OUTPUT.csv]" >&2
  exit 1
}

python3 - "$COUNT_LIST" "$OUTPUT" <<'PY'
import csv
import sys
from pathlib import Path

files = [Path(x.strip()) for x in Path(sys.argv[1]).read_text().splitlines() if x.strip()]
output = Path(sys.argv[2])

gene_ids = None
sample_counts = []
sample_names = []

for f in files:
    rows = []
    with f.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            rows.append(line.rstrip("\n").split("\t"))

    body = rows[1:]
    genes = [r[0] for r in body]
    counts = [r[-1] for r in body]

    if gene_ids is None:
        gene_ids = genes
    elif genes != gene_ids:
        raise SystemExit(f"Gene order differs in {f}")

    name = f.name
    if name.endswith("_gene_counts.txt"):
        name = name[:-len("_gene_counts.txt")]

    sample_names.append(name)
    sample_counts.append(counts)

with output.open("w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["gene_id", *sample_names])
    for i, gene in enumerate(gene_ids):
        writer.writerow([gene, *[x[i] for x in sample_counts]])

print(f"Wrote {output}")
PY
