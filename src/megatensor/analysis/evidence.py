"""Composite site evidence and cross-PXD intensity panels."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import load_site_index, pride_site_table
from megatensor.analysis.concordance import intensity_by_site_pxd
from megatensor.analysis.replication import replication_tables
from megatensor.viz.completeness import is_human_uniprot


def _valid_site_keys(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(
        pl.col("site_key").is_not_null(),
        pl.col("site_key").str.contains(r"^[^:]+:\d+:[ST]$"),
    )


def site_evidence_ladder(*, top_n: int | None = None) -> pl.DataFrame:
    """Rank sites by canon depth, PRIDE replication, and SET support."""
    idx = _valid_site_keys(load_site_index())
    pride = pride_site_table().select("site_key", "n_pxds", "pxds")
    ladder = (
        idx.join(pride, on="site_key", how="left")
        .with_columns(pl.col("n_pxds").fill_null(0))
        .with_columns(pl.col("site_key").str.split(":").list.get(0).alias("protein_acc"))
        .filter(pl.col("protein_acc").map_elements(is_human_uniprot, return_dtype=pl.Boolean))
        .with_columns(
            (
                pl.col("n_layers").cast(pl.Float64) * 10.0
                + pl.col("n_pxds").cast(pl.Float64) * 15.0
                + (pl.col("set_hits").cast(pl.Float64) + 1.0).log() * 4.0
            ).alias("evidence_score"),
        )
        .sort(["evidence_score", "n_pxds", "set_hits"], descending=True)
    )
    if top_n is not None:
        ladder = ladder.head(top_n)
    return ladder


def triangulated_intensity_matrix(*, top_n: int = 25) -> pl.DataFrame:
    """Top triangulated sites × PXD mean log₁₀ intensity (wide)."""
    tri = replication_tables()["triangulated_sites"].head(top_n)
    if tri.is_empty():
        return pl.DataFrame()

    sites = tri["site_key"].to_list()
    wide = (
        intensity_by_site_pxd()
        .filter(pl.col("site_key").is_in(sites))
        .with_columns((pl.col("intensity_mean") + 1.0).log10().alias("log10_int"))
    )
    pivot = wide.pivot(on="dataset_id", index="site_key", values="log10_int", aggregate_function="first")
    return tri.select("site_key").join(pivot, on="site_key", how="left")
