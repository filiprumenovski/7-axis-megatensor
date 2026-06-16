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

# Phase 1: PRIDE metadata snapshot + glyco discovery (stops at checkpoint)
pride-discover:
    .venv/bin/megatensor pride-discover

# Phase 2: download chosen PRIDE result tables (after Filip approves PXDs)
pride-ingest:
    .venv/bin/megatensor pride-ingest

# Phase 3: assemble Megatensor view
assemble:
    .venv/bin/megatensor assemble

# Phase 4: enrichment + figures + exports
enrich:
    .venv/bin/megatensor enrich

figures:
    .venv/bin/megatensor figures

export:
    .venv/bin/megatensor export

report:
    .venv/bin/megatensor report

# Full autonomous loop (stops at checkpoints per §12)
all: setup download canon
