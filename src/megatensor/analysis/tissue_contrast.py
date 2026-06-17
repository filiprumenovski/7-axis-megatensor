"""Tissue-resolved contrasts within single PRIDE deposits."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import load_pride_obs, site_key_expr

BAP1KO_PXD = "PXD035902"


def bap1ko_brain_liver_pairs(*, pxd: str = BAP1KO_PXD) -> pl.DataFrame:
    """Per-site mean intensity in brain vs liver (PXD035902 glycomics)."""
    obs = load_pride_obs().filter(
        pl.col("dataset_id") == pxd,
        pl.col("metric_name") == "intensity",
        pl.col("cond_tissue").is_in(["brain", "liver"]),
    )
    if obs.is_empty():
        return pl.DataFrame()

    agg = (
        obs.with_columns(site_key_expr("protein_id_raw", "residue_pos_raw").alias("site_key"))
        .group_by(["site_key", "cond_tissue"])
        .agg(pl.col("metric_value").mean().alias("intensity_mean"), pl.len().alias("n_obs"))
    )
    wide = agg.pivot(on="cond_tissue", index="site_key", values="intensity_mean")
    if "brain" not in wide.columns or "liver" not in wide.columns:
        return pl.DataFrame()

    return (
        wide.filter(pl.col("brain").is_not_null(), pl.col("liver").is_not_null(), pl.col("brain") > 0, pl.col("liver") > 0)
        .with_columns(
            (pl.col("brain") / pl.col("liver")).alias("brain_liver_ratio"),
            (pl.col("brain").log10() + 1.0).alias("log10_brain"),
            (pl.col("liver").log10() + 1.0).alias("log10_liver"),
            pl.col("site_key").str.split(":").list.get(0).alias("protein_acc"),
        )
        .sort("brain_liver_ratio", descending=True)
    )


def bap1ko_brain_liver_summary(pairs: pl.DataFrame) -> dict:
    if pairs.is_empty():
        return {}
    ratios = pairs["brain_liver_ratio"].to_numpy()
    import numpy as np

    top = pairs.row(0, named=True)
    return {
        "pxd": BAP1KO_PXD,
        "n_shared_sites": pairs.height,
        "median_brain_liver_ratio": round(float(np.median(ratios)), 2),
        "pct_brain_enriched_2x": round(100 * float((ratios >= 2).mean()), 1),
        "top_site": top["site_key"],
        "top_brain_liver_ratio": round(float(top["brain_liver_ratio"]), 1),
    }
