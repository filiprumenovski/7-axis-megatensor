"""§11.5 ML-ready tensor exports (no modeling)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl
import structlog

from megatensor.paths import EXPORTS
from megatensor.store import PRIDE_STORE, UNION_STORE

log = structlog.get_logger()


def export_site_condition_matrix(
    *,
    metric_name: str = "intensity",
    dataset_id: str | None = None,
) -> tuple[pl.DataFrame, Path, Path]:
    """Dense site x condition matrix for one metric from PRIDE layer."""
    obs = pl.read_parquet(PRIDE_STORE.staging / "observations.parquet")
    obs = obs.filter(pl.col("metric_name") == metric_name)
    if dataset_id:
        obs = obs.filter(pl.col("dataset_id") == dataset_id)
    if obs.is_empty():
        raise ValueError(f"no rows for metric_name={metric_name!r}")

    obs = obs.with_columns(
        pl.concat_str(
            [
                pl.col("protein_id_raw"),
                pl.col("residue_pos_raw").cast(pl.Utf8),
                pl.col("residue_aa"),
            ],
            separator=":",
        ).alias("site_key"),
        pl.concat_str(
            [
                pl.col("dataset_id"),
                pl.col("cond_treatment").fill_null("na"),
                pl.col("cond_tissue").fill_null("na"),
            ],
            separator="|",
        ).alias("condition_key"),
    )

    pivot = (
        obs.group_by(["site_key", "condition_key"])
        .agg(pl.col("metric_value").mean().alias("value"))
        .pivot(on="condition_key", index="site_key", values="value", aggregate_function="first")
        .fill_null(0.0)
    )

    EXPORTS.mkdir(parents=True, exist_ok=True)
    parquet_path = EXPORTS / "site_x_condition.parquet"
    pivot.write_parquet(parquet_path)

    cond_cols = [c for c in pivot.columns if c != "site_key"]
    arr = pivot.select(cond_cols).to_numpy()
    npy_path = EXPORTS / "site_x_condition.npy"
    np.save(npy_path, arr)
    meta = EXPORTS / "site_x_condition_meta.json"
    meta.write_text(
        json.dumps(
            {"site_keys": pivot["site_key"].to_list(), "condition_keys": cond_cols, "metric_name": metric_name},
            indent=2,
        )
    )
    return pivot, parquet_path, npy_path


def export_site_feature_matrix() -> tuple[pl.DataFrame, Path, Path]:
    """Enriched site x feature matrix from union enrichment."""
    feat_path = UNION_STORE.enrichment / "site_features.parquet"
    if not feat_path.is_file():
        raise FileNotFoundError("run: megatensor enrich")

    enriched = pl.read_parquet(feat_path)
    enriched = enriched.with_columns(
        pl.concat_str(
            [pl.col("protein_acc"), pl.col("residue_pos").cast(pl.Utf8), pl.col("residue_aa")],
            separator=":",
        ).alias("site_key"),
        pl.col("region_type").fill_null("none").alias("region_type_filled"),
    )

    one_hot = enriched.to_dummies(columns=["region_type_filled"], separator="=")
    feature_cols = [c for c in one_hot.columns if c not in ("identity_id", "protein_acc", "residue_pos", "residue_aa", "site_key", "domain_name", "seq_window", "gene_symbol", "seq_match")]
    base = ["disorder_score"]
    use_cols = [c for c in base if c in one_hot.columns] + [c for c in one_hot.columns if c.startswith("region_type_filled=")]

    mat = one_hot.select(["site_key", *use_cols]).fill_null(0.0)
    EXPORTS.mkdir(parents=True, exist_ok=True)
    parquet_path = EXPORTS / "site_x_features.parquet"
    mat.write_parquet(parquet_path)

    arr = mat.drop("site_key").to_numpy()
    npy_path = EXPORTS / "site_x_features.npy"
    np.save(npy_path, arr)
    (EXPORTS / "site_x_features_meta.json").write_text(
        json.dumps({"site_keys": mat["site_key"].to_list(), "feature_columns": use_cols}, indent=2)
    )
    return mat, parquet_path, npy_path


def run_export() -> dict:
    site_cond, p1, n1 = export_site_condition_matrix(metric_name="intensity")
    site_feat, p2, n2 = export_site_feature_matrix()

    doc = f"""# Megatensor exports

ML-ready matrices (no modeling in this repo).

## site_x_condition

- Parquet: `{p1}`
- NumPy: `{n1}`
- Shape: {site_cond.height} sites x {site_cond.width - 1} conditions
- Source: PRIDE observation layer, `metric_name=intensity`

## site_x_features

- Parquet: `{p2}`
- NumPy: `{n2}`
- Shape: {site_feat.height} sites x {site_feat.width - 1} features
- Source: union enrichment (`site_features.parquet`)

Join on `site_key` = `protein_acc:residue_pos:residue_aa`.
"""
    EXPORTS.mkdir(parents=True, exist_ok=True)
    (EXPORTS / "EXPORTS.md").write_text(doc)

    summary = {
        "site_x_condition": {"sites": site_cond.height, "conditions": site_cond.width - 1},
        "site_x_features": {"sites": site_feat.height, "features": site_feat.width - 1},
    }
    log.info("export_complete", **summary)
    print(json.dumps(summary, indent=2))
    return summary
