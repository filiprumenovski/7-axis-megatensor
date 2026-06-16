"""§9 SET assembly and Megatensor append."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from megatensor.hash_utils import hash_u64
from megatensor.paths import METRICS, REGISTRY, SETS


def assign_set_uid(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        pl.struct(
            [
                "identity_id",
                "ptm_id",
                "condition_id",
                "acquisition_id",
                "instrument_id",
                "provenance_id",
            ]
        )
        .map_elements(
            lambda r: hash_u64(
                r["identity_id"],
                r["ptm_id"],
                r["condition_id"],
                r["acquisition_id"],
                r["instrument_id"],
                r["provenance_id"],
            ),
            return_dtype=pl.UInt64,
        )
        .alias("set_uid")
    )


def write_partition(
    coords: pl.DataFrame,
    metrics: pl.DataFrame,
    dataset_id: str,
) -> None:
    coord_path = SETS / f"dataset_id={dataset_id}"
    metric_path = METRICS / f"dataset_id={dataset_id}"
    coord_path.mkdir(parents=True, exist_ok=True)
    metric_path.mkdir(parents=True, exist_ok=True)

    coord_cols = [
        "set_uid",
        "identity_id",
        "ptm_id",
        "condition_id",
        "acquisition_id",
        "instrument_id",
        "provenance_id",
        "dataset_id",
        "loc_score",
        "loc_method",
        "loc_is_ambiguous",
    ]
    coords.select(coord_cols).unique(subset=["set_uid"]).write_parquet(coord_path / "part.parquet")
    metrics.select(
        "set_uid",
        "metric_name",
        "metric_value",
        "metric_norm_state",
        "metric_unit",
        "qc_flags",
        "dataset_id",
    ).write_parquet(metric_path / "part.parquet")


def write_registry(identity_dim: pl.DataFrame, dims: dict[str, pl.DataFrame]) -> None:
    REGISTRY.mkdir(parents=True, exist_ok=True)
    identity_dim.write_parquet(REGISTRY / "identity_dim.parquet")
    for name, table in dims.items():
        table.write_parquet(REGISTRY / f"{name}.parquet")
