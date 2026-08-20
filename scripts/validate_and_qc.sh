#!/usr/bin/env bash
# validate_and_qc.sh
# Step 3 (confirm FASTA files are correct) + Step 4 (assembly-quality checks)
# from the project instructions. Requires: seqkit (conda install -c bioconda seqkit)
#
# Usage:
#   bash scripts/validate_and_qc.sh

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GENOME_DIR="${ROOT_DIR}/data/genomes"
META="${ROOT_DIR}/data/metadata/samplesheet.tsv"
OUTDIR="${ROOT_DIR}/results/qc"
mkdir -p "${OUTDIR}"

echo "== Step 3: basic FASTA sanity checks =="
fail=0
for fq in "${GENOME_DIR}"/*.fna.gz; do
  [[ -e "${fq}" ]] || { echo "No genome files found in ${GENOME_DIR}. Run download_genomes.sh first." >&2; exit 1; }
  name=$(basename "${fq}")

  # not empty
  size=$(stat -c%s "${fq}" 2>/dev/null || stat -f%z "${fq}")
  if [[ "${size}" -lt 100 ]]; then
    echo "  FAIL: ${name} is empty or too small"; fail=1; continue
  fi

  # must be nucleotide (ACGTN), not protein
  non_nt=$(zcat "${fq}" | grep -v '^>' | tr -d 'ACGTNacgtn\n' | wc -c)
  if [[ "${non_nt}" -gt 0 ]]; then
    echo "  WARN: ${name} contains ${non_nt} non-ACGTN characters (check it's nucleotide, not protein)"
  fi

  # name matches metadata sample_id
  sample_id="${name%.fna.gz}"
  if ! grep -q "^${sample_id}"$'\t' "${META}"; then
    echo "  FAIL: ${sample_id} not listed in samplesheet.tsv"; fail=1; continue
  fi

  echo "  OK: ${name}"
done

echo ""
echo "== Step 4: assembly-quality stats (seqkit) =="
if command -v seqkit >/dev/null 2>&1; then
  seqkit stats -a "${GENOME_DIR}"/*.fna.gz | tee "${OUTDIR}/seqkit_stats.tsv"
else
  echo "  seqkit not found on PATH."
  echo "  Install with: conda install -c bioconda -c conda-forge seqkit"
  echo "  (or: mamba install -c bioconda seqkit)"
fi

echo ""
echo "== Duplicate check (by md5 of decompressed content) =="
md5_report="${OUTDIR}/md5_check.tsv"
: > "${md5_report}"
for fq in "${GENOME_DIR}"/*.fna.gz; do
  sum=$(zcat "${fq}" | md5sum | cut -d' ' -f1)
  echo -e "$(basename "${fq}")\t${sum}" >> "${md5_report}"
done
dupes=$(cut -f2 "${md5_report}" | sort | uniq -d)
if [[ -n "${dupes}" ]]; then
  echo "  WARNING: duplicate genome content detected (identical md5):"
  echo "${dupes}"
else
  echo "  OK: no duplicate genome content detected."
fi

if [[ "${fail}" -eq 1 ]]; then
  echo ""
  echo "One or more checks FAILED. Fix before proceeding to Step 5." >&2
  exit 1
fi

echo ""
echo "All Step 3 checks passed. Review ${OUTDIR}/seqkit_stats.tsv before Step 5 (species confirmation)."
