"""Chemoproteomics and genetic perturbation summaries."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import load_pride_obs


def chemoproteomics_site_matrix(pxd: str = "PXD063995") -> pl.DataFrame:
    """Site × probe matrix (spectral support or intensity per probe)."""
    obs = load_pride_obs().filter(
        pl.col("dataset_id") == pxd,
        pl.col("cond_treatment").is_not_null(),
    )
    if obs.is_empty():
        return pl.DataFrame()

    metric = "intensity" if obs.filter(pl.col("metric_name") == "intensity").height else "spectral_count"
    obs = obs.filter(pl.col("metric_name") == metric)

    return (
        obs.with_columns(
            pl.concat_str(
                [pl.col("protein_id_raw"), pl.col("residue_pos_raw").cast(pl.Utf8), pl.col("residue_aa")],
                separator=":",
            ).alias("site_key")
        )
        .group_by(["site_key", "cond_treatment"])
        .agg(pl.col("metric_value").sum().alias("value"))
        .pivot(on="cond_treatment", index="site_key", values="value")
        .fill_null(0.0)
    )


def bap1ko_tissue_sites(pxd: str = "PXD035902") -> pl.DataFrame:
    obs = load_pride_obs().filter(pl.col("dataset_id") == pxd)
    return (
        obs.with_columns(
            pl.concat_str(
                [pl.col("protein_id_raw"), pl.col("residue_pos_raw").cast(pl.Utf8), pl.col("residue_aa")],
                separator=":",
            ).alias("site_key")
        )
        .group_by(["site_key", "cond_tissue", "cond_treatment"])
        .agg(
            pl.col("metric_value").filter(pl.col("metric_name") == "intensity").mean().alias("intensity"),
            pl.len().alias("n_obs"),
        )
        .sort("n_obs", descending=True)
    )
