#!/usr/bin/env bash

set -euo pipefail

COMPRESS_ONLY=0
if [[ "${1:-}" == "--compress-only" ]]; then
    COMPRESS_ONLY=1
    shift
fi

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
XLSX="${BASE_DIR}/data/WHO-UCN-TB-2023.xlsx"
PREFIX="${1:-${BASE_DIR}/data/$(date +%Y-%m-%d)}"

TSV="${PREFIX}_annotation_table.tsv"
HEADER="${PREFIX}_annotations_header.txt"

for tool in bgzip tabix python3; do
    command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: ${tool} not on PATH" >&2; exit 1; }
done

if [[ "${COMPRESS_ONLY}" -eq 1 ]]; then
    if [[ ! -f "${TSV}" ]]; then
        echo "ERROR: --compress-only needs an existing ${TSV}" >&2
        exit 1
    fi
    echo "==> Skipping spreadsheet parse; compressing existing ${TSV}"
else
    if [[ ! -f "${XLSX}" ]]; then
        echo "ERROR: WHO catalogue not found at ${XLSX}" >&2
        exit 1
    fi

    COORD_VCF="$(ls -1 "${BASE_DIR}"/data/Genomic_coordinates_*.vcf.gz 2>/dev/null | head -1 || true)"
    COORD_ARG=()
    if [[ -n "${COORD_VCF}" ]]; then
        echo "==> Coordinate cross-check against $(basename "${COORD_VCF}")"
        COORD_ARG=(--coordinates_vcf "${COORD_VCF}")
    else
        echo "==> No Genomic_coordinates_*.vcf.gz found; skipping coordinate cross-check"
    fi

    echo "==> Building annotation table from ${XLSX}"
    python3 "${BASE_DIR}/tools/annotation/prepare_who_annotation.py" \
        --input  "${XLSX}" \
        --output "${TSV}" \
        --header "${HEADER}" \
        "${COORD_ARG[@]}"
fi

echo "==> Compressing and indexing"
bgzip -f "${TSV}"
tabix -f -s1 -b2 -e2 "${TSV}.gz"

echo "==> Verifying the parallel-array invariant"
python3 - "${TSV}.gz" <<'PYCHECK'
import gzip, sys
from collections import Counter

path = sys.argv[1]
allele_rows = multi = bad = placeholder = unsorted = 0
example = None
prev_key = None
grades = Counter()

with gzip.open(path, 'rt') as handle:
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 10:
            continue
        allele_rows += 1

        key = (fields[0], int(fields[1]))
        if prev_key is not None and key < prev_key:
            unsorted += 1
        prev_key = key

        if fields[5] == 'unknown':
            continue

        drugs  = fields[5].split(',')
        gr     = fields[7].split(',')
        ids    = fields[8].split(',')

        if len(drugs) > 1:
            multi += 1
        if not (len(drugs) == len(gr) == len(ids)):
            bad += 1
            if example is None:
                example = line.rstrip('\n')
        if any(g.strip().lower() == 'unknown' for g in gr):
            placeholder += 1
        grades.update(g.strip() for g in gr)

print(f"    allele rows           : {allele_rows}")
print(f"    multi-drug alleles    : {multi}")
print(f"    misaligned arrays     : {bad}")
print(f"    'unknown' in a grade  : {placeholder}")
print(f"    out-of-order rows     : {unsorted}")
print("    grade distribution (variant-drug pairs):")
for grade, count in grades.most_common():
    print(f"      {grade:<26} {count}")

if bad or placeholder or unsorted:
    if example:
        print(f"    first misaligned row  : {example}")
    print("    FAILED: the table is not safe for per-drug interpretation.")
    sys.exit(1)
PYCHECK

echo
echo "Done."
echo "  table  : ${TSV}.gz"
echo "  index  : ${TSV}.gz.tbi"
echo "  header : ${HEADER}"
echo
echo "Now set in nextflow.config:"
echo "  annotation_table  = \"\$baseDir/data/$(basename "${PREFIX}")_annotation_table.tsv.gz\""
echo "  annotation_header = \"\$baseDir/data/$(basename "${PREFIX}")_annotations_header.txt\""
