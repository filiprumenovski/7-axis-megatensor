"""Phase 2b: parse PRIDE downloads -> observation rows + QC."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import structlog

from megatensor.ingest.pride_router import parse_downloads
from megatensor.paths import FIGURES, PRIDE_DOWNLOADS
from megatensor.store import PRIDE_STORE

log = structlog.get_logger()


def _qc_summary(obs: pl.DataFrame, qc: pl.DataFrame) -> dict:
    by_pxd = (
        obs.group_by("dataset_id")
        .agg(
            pl.len().alias("observation_rows"),
            pl.struct(["protein_id_raw", "residue_pos_raw", "residue_aa"]).n_unique().alias("unique_sites"),
            pl.col("protein_id_raw").n_unique().alias("unique_proteins"),
            pl.col("metric_name").unique().alias("metric_kinds"),
            pl.col("metric_unit").drop_nulls().unique().alias("metric_units"),
            pl.col("source_engine").unique().alias("engines"),
            pl.col("loc_score").null_count().alias("loc_score_nulls"),
            pl.col("cond_tissue").is_not_null().mean().alias("pct_tissue"),
            pl.col("cond_treatment").is_not_null().mean().alias("pct_treatment"),
            pl.col("acq_collision").is_not_null().mean().alias("pct_collision"),
            pl.col("inst_ms_cv").is_not_null().mean().alias("pct_inst_cv"),
        )
        if obs.height
        else pl.DataFrame()
    )
    by_engine_metric = (
        obs.group_by(["source_engine", "metric_name", "metric_unit"])
        .agg(pl.len().alias("rows"), pl.col("metric_value").median().alias("median_value"))
        .sort("rows", descending=True)
        if obs.height
        else pl.DataFrame()
    )
    return {
        "files_parsed": int(qc.height),
        "files_ok": int(qc.filter(pl.col("status") == "ok").height),
        "files_empty": int(qc.filter(pl.col("status") == "empty").height),
        "files_missing": int(qc.filter(pl.col("status") == "missing").height),
        "files_error": int(qc.filter(pl.col("status").str.starts_with("error")).height),
        "total_observation_rows": int(obs.height),
        "total_unique_sites": int(
            obs.select(["protein_id_raw", "residue_pos_raw", "residue_aa"]).unique().height
        )
        if obs.height
        else 0,
        "by_pxd": by_pxd.to_dicts() if by_pxd.height else [],
        "by_engine_metric": by_engine_metric.to_dicts() if by_engine_metric.height else [],
        "engines": obs.group_by("source_engine").len().to_dicts() if obs.height else [],
    }


def run_pride_ingest(download_root: Path | None = None) -> dict:
    root = download_root or PRIDE_DOWNLOADS
    if not root.exists():
        raise FileNotFoundError(f"No downloads at {root} — run: just pride-download")

    obs, qc = parse_downloads(root)
    obs.write_parquet(PRIDE_STORE.staging / "observations.parquet")
    qc.write_csv(PRIDE_STORE.staging / "parse_qc.csv")

    FIGURES.mkdir(parents=True, exist_ok=True)
    summary = _qc_summary(obs, qc)
    (FIGURES / "pride_ingest_qc.json").write_text(json.dumps(summary, indent=2))
    if obs.height:
        obs.group_by(["source_engine", "metric_name", "metric_unit", "dataset_id"]).len().sort(
            "len", descending=True
        ).write_csv(FIGURES / "pride_engine_metrics.csv")
    qc.write_csv(FIGURES / "pride_ingest_files.csv")

    log.info("pride_ingest_complete", **{k: v for k, v in summary.items() if k != "by_pxd"})
    print(json.dumps(summary, indent=2))
    return summary
