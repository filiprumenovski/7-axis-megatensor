#!/usr/bin/env bash
# Unpack a pride-ingest parquet snapshot tarball into data/pride/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARCHIVE="${1:-${ROOT}/pride_snapshot_parquets_2026-04-07.tar.gz}"
DEST="${ROOT}/data/pride"

if [[ ! -f "${ARCHIVE}" ]]; then
  echo "archive not found: ${ARCHIVE}" >&2
  exit 1
fi

mkdir -p "${DEST}"
echo "extracting ${ARCHIVE} -> ${DEST}"
tar -xzf "${ARCHIVE}" -C "${DEST}" --strip-components=1 \
  --exclude='._*' --exclude='*.DS_Store' 2>/dev/null || \
tar -xzf "${ARCHIVE}" -C "${ROOT}/data" --exclude='._*' --exclude='*.DS_Store'

# tarball may unpack as data/bronze or bronze directly
if [[ -d "${ROOT}/data/bronze" && ! -d "${DEST}/bronze" ]]; then
  rsync -a "${ROOT}/data/bronze" "${ROOT}/data/silver" "${DEST}/"
fi

echo "snapshot ready:"
du -sh "${DEST}"/* 2>/dev/null || true
