#!/usr/bin/env bash
# Bulk-download canonical O-GlcNAc reference datasets (no API calls).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${ROOT}/data/canon"
mkdir -p "${DEST}"

MCW_URL='https://www.oglcnac.mcw.edu/download/?oglcnac_organisms=All%20species&oglcnac_format=csv&download_oglcnac=download_oglcnac'
ATLAS_I='https://oglcnac.org/static/dataset/Atlas%205.0_unambiguous%20sites_20251208.csv'
ATLAS_II='https://oglcnac.org/static/dataset/Atlas%205.0_ambiguous%20sites_20251208.csv'

download() {
  local url="$1" out="$2"
  if [[ -f "${out}" ]]; then
    echo "exists: ${out}"
    return 0
  fi
  echo "downloading: ${out}"
  curl -fL --retry 3 --retry-delay 5 -o "${out}.part" "${url}"
  mv "${out}.part" "${out}"
}

# MCW site requires -k on some trust chains
download_mcw() {
  local out="${DEST}/oglcnacdb_all_species.csv"
  if [[ -f "${out}" ]]; then
    echo "exists: ${out}"
    return 0
  fi
  echo "downloading: ${out}"
  curl -fkL --retry 3 --retry-delay 5 -o "${out}.part" "${MCW_URL}"
  mv "${out}.part" "${out}"
}

download_mcw
download "${ATLAS_I}"  "${DEST}/atlas_dataset_I_unambiguous.csv"
download "${ATLAS_II}" "${DEST}/atlas_dataset_II_ambiguous.csv"

echo "canon bulk files in ${DEST}:"
ls -lh "${DEST}"
