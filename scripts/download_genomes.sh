#!/usr/bin/env bash
# download_genomes.sh
# Downloads the 16 selected Bacteroides vulgatus assemblies from CNGB (CNP0000126),
# renames them to the Bvulgatus_NNN.fna.gz convention, and verifies each file.
#
# Usage:
#   bash scripts/download_genomes.sh
#
# Requires: wget (or curl), gzip. Reads data/metadata/samplesheet.tsv.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SAMPLESHEET="${ROOT_DIR}/data/metadata/samplesheet.tsv"
OUTDIR="${ROOT_DIR}/data/genomes"

mkdir -p "${OUTDIR}"

if [[ ! -f "${SAMPLESHEET}" ]]; then
  echo "ERROR: samplesheet not found at ${SAMPLESHEET}" >&2
  exit 1
fi

echo "== Downloading genomes listed in ${SAMPLESHEET} =="

# Skip header, read tab-separated columns: sample_id species accession source original_name fasta_file_name ftp_url ...
tail -n +2 "${SAMPLESHEET}" | while IFS=$'\t' read -r sample_id species accession source original_name fasta_file_name ftp_url contigs n50 total_size; do
  dest="${OUTDIR}/${sample_id}.fna.gz"

  if [[ -s "${dest}" ]]; then
    echo "  [skip] ${sample_id} already present"
    continue
  fi

  echo "  [get]  ${sample_id}  <-  ${ftp_url}"
  # -q quiet, retry on transient FTP failures
  if ! wget -q --tries=3 --timeout=60 -O "${dest}.tmp" "${ftp_url}"; then
    echo "  ERROR: failed to download ${ftp_url}" >&2
    rm -f "${dest}.tmp"
    continue
  fi

  # sanity check: file must be a valid gzip and non-empty
  if ! gzip -t "${dest}.tmp" 2>/dev/null; then
    echo "  ERROR: ${sample_id} downloaded file is not valid gzip" >&2
    rm -f "${dest}.tmp"
    continue
  fi

  mv "${dest}.tmp" "${dest}"
  echo "  [ok]   ${sample_id} -> $(basename "${dest}")"
done

echo "== Download complete. Files in ${OUTDIR}: =="
ls -la "${OUTDIR}"/*.fna.gz 2>/dev/null || echo "  (no files found - check errors above)"
