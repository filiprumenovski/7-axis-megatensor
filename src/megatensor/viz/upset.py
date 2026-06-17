"""Figure B: UpSet-ready site membership exports."""

from __future__ import annotations

import polars as pl

from megatensor.store import CANON_STORE, PRIDE_STORE
from megatensor.viz.completeness import canon_overlap


def _site_keys(store, layer: str) -> pl.DataFrame:
    identity = pl.read_parquet(store.registry / "identity_dim.parquet")
    coords = pl.read_parquet(str(store.sets / "dataset_id=*/*.parquet"), hive_partitioning=True)
    joined = coords.join(identity, on="identity_id", how="left").filter(~pl.col("protein_level_only"))
    return joined.with_columns(
        pl.lit(layer).alias("layer"),
        pl.concat_str(
            [pl.col("protein_acc"), pl.col("residue_pos").cast(pl.Utf8), pl.col("residue_aa")],
            separator=":",
        ).alias("site_key"),
    )


def export_upset_tables() -> dict:
    canon = _site_keys(CANON_STORE, "canon")
    pride = _site_keys(PRIDE_STORE, "pride")

    canon_membership = (
        canon.select("site_key", "dataset_id")
        .unique()
        .pivot(on="dataset_id", index="site_key", values="dataset_id", aggregate_function="len")
        .fill_null(0)
    )
    pride_membership = (
        pride.select("site_key", "dataset_id")
        .unique()
        .pivot(on="dataset_id", index="site_key", values="dataset_id", aggregate_function="len")
        .fill_null(0)
    )

    return {
        "canon_overlap_pairs": canon_overlap(canon),
        "canon_membership": canon_membership,
        "pride_membership": pride_membership,
    }
