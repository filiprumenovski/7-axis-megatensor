"""BAP1KO × tissue site burden (PXD035902)."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.conditions import bap1ko_tissue_sites


def bap1ko_tissue_summary(pxd: str = "PXD035902") -> pl.DataFrame:
    raw = bap1ko_tissue_sites(pxd=pxd)
    if raw.is_empty():
        return pl.DataFrame()
    return (
        raw.filter(pl.col("cond_tissue").is_not_null())
        .group_by("cond_tissue")
        .agg(
            pl.col("site_key").n_unique().alias("n_sites"),
            pl.col("n_obs").sum().alias("n_observations"),
        )
        .sort("n_sites", descending=True)
    )


def bap1ko_top_sites(*, top_n: int = 30, pxd: str = "PXD035902") -> pl.DataFrame:
    raw = bap1ko_tissue_sites(pxd=pxd)
    if raw.is_empty():
        return pl.DataFrame()
    return (
        raw.group_by("site_key", "cond_tissue")
        .agg(pl.col("n_obs").sum().alias("n_obs"))
        .sort("n_obs", descending=True)
        .head(top_n)
    )
