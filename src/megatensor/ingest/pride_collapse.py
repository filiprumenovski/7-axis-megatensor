"""Collapse redundant PRIDE observation rows before tensorization."""

from __future__ import annotations

import polars as pl

from megatensor.ingest.pride_common import finalize_obs
from megatensor.schema import OBSERVATION_COLUMNS

METRIC_FIELDS = ("metric_name", "metric_value", "metric_norm_state", "metric_unit")
GROUP_COLS = [
    c
    for c in OBSERVATION_COLUMNS
    if c
    not in (
        *METRIC_FIELDS,
        "qc_flags",
        "loc_score",
        "loc_method",
        "loc_is_ambiguous",
    )
]

# Engines where cond_replicate holds spectrum / PSM id (not a true replicate axis).
_SPECTRUM_REPLICATE_ENGINES = frozenset({"pd", "mztab"})


def _collapse_metrics(df: pl.DataFrame, group_cols: list[str]) -> pl.DataFrame:
    if df.is_empty():
        return df

    def _agg(name: str, how: str) -> pl.DataFrame:
        sub = df.filter(pl.col("metric_name") == name)
        if sub.is_empty():
            return sub
        return sub.group_by(group_cols + ["metric_name"]).agg(
            getattr(pl.col("metric_value"), how)().alias("metric_value"),
            pl.col("metric_norm_state").first(),
            pl.col("metric_unit").first(),
        )

    parts = [_agg("qvalue", "min"), _agg("spectral_count", "sum")]
    rest = df.filter(~pl.col("metric_name").is_in(["qvalue", "spectral_count"]))
    if rest.height:
        parts.append(
            rest.group_by(group_cols + ["metric_name"]).agg(
                pl.col("metric_value").max(),
                pl.col("metric_norm_state").first(),
                pl.col("metric_unit").first(),
            )
        )
    return pl.concat([p for p in parts if p.height > 0], how="vertical_relaxed")


def _collapse_group(df: pl.DataFrame, group_cols: list[str]) -> pl.DataFrame:
    if df.is_empty():
        return df
    metrics = _collapse_metrics(df, group_cols)
    meta = df.group_by(group_cols).agg(
        pl.col("loc_score").max(),
        pl.col("loc_method").first(),
        pl.col("loc_is_ambiguous").any(),
        pl.col("qc_flags").first(),
    )
    return finalize_obs(metrics.join(meta, on=group_cols, how="left"))


def collapse_pride_observations(df: pl.DataFrame) -> pl.DataFrame:
    """Roll up PSM/spectrum fan-out and deduplicate metrics per site event."""
    if df.is_empty():
        return df
    df = df.with_columns(
        pl.when(pl.col("source_engine").is_in(list(_SPECTRUM_REPLICATE_ENGINES)))
        .then(None)
        .otherwise(pl.col("cond_replicate"))
        .alias("cond_replicate")
    )
    return _collapse_group(df, GROUP_COLS)
