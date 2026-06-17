"""Cross-study site replication and canon triangulation."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import load_site_index, pride_site_table
from megatensor.viz.completeness import is_human_uniprot


def replication_tables() -> dict[str, pl.DataFrame]:
    idx = load_site_index()
    pride = pride_site_table()

    multi_pxd = pride.filter(pl.col("n_pxds") >= 2).sort("n_pxds", descending=True)
    triangulated = (
        idx.filter(pl.col("n_layers") >= 2)
        .join(pride.select("site_key", "n_pxds", "pxds"), on="site_key", how="left")
        .filter(pl.col("n_pxds") >= 2)
        .sort("n_pxds", descending=True)
    )

    human_multi = multi_pxd.filter(pl.col("protein_acc").map_elements(is_human_uniprot, return_dtype=pl.Boolean))

    protein_hubs = (
        triangulated.with_columns(pl.col("site_key").str.split(":").list.get(0).alias("protein_acc"))
        .group_by("protein_acc")
        .agg(
            pl.len().alias("n_triangulated_sites"),
            pl.col("n_pxds").max().alias("max_pxds"),
            pl.col("site_key").head(5).alias("example_sites"),
        )
        .sort("n_triangulated_sites", descending=True)
    )

    overlap_tiers = pl.DataFrame(
        {
            "tier": [
                "pride_unique",
                "pride_multi_pxd",
                "canon_and_pride",
                "canon_and_pride_multi_pxd",
            ],
            "n_sites": [
                idx.filter(pl.col("n_layers") == 1, pl.col("layers").list.contains("pride")).height,
                multi_pxd.height,
                idx.filter(pl.col("n_layers") >= 2).height,
                triangulated.height,
            ],
        }
    )

    return {
        "multi_pxd_sites": multi_pxd,
        "triangulated_sites": triangulated,
        "human_multi_pxd_sites": human_multi,
        "protein_hubs": protein_hubs,
        "overlap_tiers": overlap_tiers,
    }
