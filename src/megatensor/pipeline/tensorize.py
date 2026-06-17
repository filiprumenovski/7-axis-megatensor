"""Shared observation -> SET tensorization for an isolated store."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.context import encode_context
from megatensor.paths import FIGURES
from megatensor.registry import resolve_identity
from megatensor.sets import assign_set_uid, write_partition, write_registry
from megatensor.store import TensorStore
from megatensor.viz.completeness import axis_completeness

log = structlog.get_logger()


def tensorize_observations(
    obs: pl.DataFrame,
    store: TensorStore,
    *,
    extra_summary: dict | None = None,
) -> tuple[dict, pl.DataFrame]:
    """Resolve identity, encode context axes, write SET partitions into one store."""
    store.root.mkdir(parents=True, exist_ok=True)
    store.staging.mkdir(parents=True, exist_ok=True)
    obs.write_parquet(store.staging / "observations.parquet")

    resolved, identity_dim = resolve_identity(obs)
    encoded, dims = encode_context(resolved)
    with_uid = assign_set_uid(encoded)
    write_registry(identity_dim, dims, store=store)

    for dataset_id in with_uid["dataset_id"].unique().to_list():
        part = with_uid.filter(pl.col("dataset_id") == dataset_id)
        write_partition(part, part, str(dataset_id), store=store)

    FIGURES.mkdir(parents=True, exist_ok=True)
    axis_completeness(with_uid).write_csv(FIGURES / f"axis_completeness_{store.layer}.csv")

    summary = {
        "layer": store.layer,
        "observation_rows": obs.height,
        "unique_sites": identity_dim.filter(~pl.col("protein_level_only")).height,
        "unique_sets": with_uid["set_uid"].n_unique(),
        "datasets": with_uid["dataset_id"].unique().sort().to_list(),
    }
    if extra_summary:
        summary.update(extra_summary)
    store.summary_path.write_text(json.dumps(summary, indent=2))
    log.info(f"{store.layer}_tensorized", **summary)
    return summary, with_uid
