"""SILAC Light/Heavy analysis (PXD039536)."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import load_pride_obs

SILAC_PXD = "PXD039536"


def silac_fold_changes(pxd: str = SILAC_PXD) -> pl.DataFrame:
    obs = load_pride_obs().filter(
        pl.col("dataset_id") == pxd,
        pl.col("metric_name") == "intensity",
        pl.col("cond_treatment").is_in(["Light", "Heavy"]),
    )
    if obs.is_empty():
        return pl.DataFrame()

    agg = obs.group_by(["protein_id_raw", "residue_pos_raw", "residue_aa"]).agg(
        pl.col("metric_value").filter(pl.col("cond_treatment") == "Heavy").mean().alias("heavy_mean"),
        pl.col("metric_value").filter(pl.col("cond_treatment") == "Light").mean().alias("light_mean"),
        pl.col("metric_value").filter(pl.col("cond_treatment") == "Heavy").len().alias("n_heavy"),
        pl.col("metric_value").filter(pl.col("cond_treatment") == "Light").len().alias("n_light"),
        pl.col("loc_score").max().alias("loc_score"),
    )
    return (
        agg.filter(pl.col("light_mean") > 0, pl.col("heavy_mean") > 0)
        .with_columns(
            pl.concat_str(
                [pl.col("protein_id_raw"), pl.col("residue_pos_raw").cast(pl.Utf8), pl.col("residue_aa")],
                separator=":",
            ).alias("site_key"),
            (pl.col("heavy_mean") / pl.col("light_mean")).alias("fc"),
        )
        .with_columns(
            pl.col("fc").log(base=2).alias("log2_fc"),
            ((pl.col("heavy_mean") * pl.col("light_mean")).sqrt().log(base=2)).alias("log2_a"),
        )
        .sort("log2_fc", descending=True)
    )


def silac_summary(fc: pl.DataFrame) -> dict:
    if fc.is_empty():
        return {}
    arr = fc["log2_fc"].to_numpy()
    return {
        "pxd": SILAC_PXD,
        "n_sites": fc.height,
        "median_log2_fc": float(fc["log2_fc"].median()),
        "median_fc": float(fc["fc"].median()),
        "pct_heavy_biased": round(100 * (fc["log2_fc"] > 0).mean(), 1),
        "top_site": fc.row(0, named=True)["site_key"],
        "top_log2_fc": round(float(fc["log2_fc"][0]), 2),
    }
