# FEHL Megatensor — one-command runners (see FEHL_MEGATENSOR_BUILD.md §12)

set shell := ["bash", "-euo", "pipefail", "-c"]

default:
    @just --list

# Bootstrap venv + editable install
setup:
    python3 -m venv .venv
    .venv/bin/pip install -U pip wheel
    .venv/bin/pip install -e ".[dev]"

# Bulk-download canonical reference files (no API)
download:
    bash scripts/download_canon.sh

# Phase 0: canon spine -> observations -> registry -> SETs
canon:
    .venv/bin/megatensor canon

# Phase 1: glyco discovery from local snapshot (default)
unpack-pride:
    bash scripts/unpack_pride_snapshot.sh

pride-discover:
    .venv/bin/megatensor pride-discover

pride-discover-live:
    .venv/bin/megatensor pride-discover --live-ingest

# Phase 2: download + parse PRIDE result tables
pride-download:
    .venv/bin/megatensor pride-download

pride-download-dry:
    .venv/bin/megatensor pride-download --dry-run

pride-ingest:
    .venv/bin/megatensor pride-ingest

pride-tensorize:
    .venv/bin/megatensor pride-tensorize

union:
    .venv/bin/megatensor union

# deprecated alias
assemble:
    .venv/bin/megatensor pride-tensorize

# Phase 4: enrichment + figures + exports
enrich:
    .venv/bin/megatensor enrich

figures:
    .venv/bin/megatensor figures

export:
    .venv/bin/megatensor export

report:
    .venv/bin/megatensor report

# Panel garnish: PNG figures + enrichment + refreshed report
panel:
    .venv/bin/megatensor panel

# Phase 5: analysis + bioRxiv preprint
analyze:
    .venv/bin/megatensor analyze

analysis-figures:
    .venv/bin/megatensor analysis-figures

biorxiv:
    .venv/bin/megatensor biorxiv

publish:
    .venv/bin/megatensor publish

# Full autonomous loop (stops at checkpoints per §12)
all: setup download canon pride-discover pride-download pride-tensorize union figures enrich export report

# Phase 4 finalize (after pride-tensorize)
finalize: union figures enrich export report
