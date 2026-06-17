"""Shared helpers for megatensor analyses."""

from __future__ import annotations

import polars as pl

from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE
from megatensor.viz.completeness import is_human_uniprot


def site_key_expr(acc: str, pos: str, aa: str = "residue_aa") -> pl.Expr:
    return pl.concat_str([pl.col(acc), pl.col(pos).cast(pl.Utf8), pl.col(aa)], separator=":")


def load_site_index() -> pl.DataFrame:
    return pl.read_parquet(UNION_STORE.staging / "site_index.parquet")


def load_pride_obs() -> pl.DataFrame:
    return pl.read_parquet(PRIDE_STORE.staging / "observations.parquet")


def load_identity(layer: str) -> pl.DataFrame:
    store = CANON_STORE if layer == "canon" else PRIDE_STORE
    return pl.read_parquet(store.registry / "identity_dim.parquet")


def pride_site_table() -> pl.DataFrame:
    """One row per PRIDE site with dataset list."""
    identity = load_identity("pride").filter(~pl.col("protein_level_only"))
    coords = pl.read_parquet(str(PRIDE_STORE.sets / "dataset_id=*/*.parquet"), hive_partitioning=True)
    return (
        coords.join(identity, on="identity_id", how="left")
        .with_columns(site_key_expr("protein_acc", "residue_pos").alias("site_key"))
        .group_by("site_key", "protein_acc", "residue_pos", "residue_aa")
        .agg(
            pl.col("dataset_id").unique().sort().alias("pxds"),
            pl.col("dataset_id").n_unique().alias("n_pxds"),
            pl.col("loc_score").max().alias("loc_score_max"),
        )
    )


def human_site_keys(keys: pl.Series | list[str]) -> list[str]:
    if isinstance(keys, pl.Series):
        keys = keys.to_list()
    return [k for k in keys if is_human_uniprot(k.split(":")[0] if ":" in k else k)]
