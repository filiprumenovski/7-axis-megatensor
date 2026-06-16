"""Phase 0 canon pipeline orchestration."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import structlog

from megatensor.context import encode_context
from megatensor.ingest.canon_atlas import AtlasDatasetIAdapter, AtlasDatasetIIAdapter
from megatensor.ingest.canon_oglcnacdb import OGlcnacDbAdapter
from megatensor.paths import CANON, FIGURES, MT
from megatensor.registry import resolve_identity
from megatensor.sets import assign_set_uid, write_partition, write_registry
from megatensor.viz.completeness import axis_completeness, canon_overlap

log = structlog.get_logger()


def run_canon() -> None:
    CANON.mkdir(parents=True, exist_ok=True)
    MT.mkdir(parents=True, exist_ok=True)

    mcw = CANON / "oglcnacdb_all_species.csv"
    atlas_i = CANON / "atlas_dataset_I_unambiguous.csv"
    atlas_ii = CANON / "atlas_dataset_II_ambiguous.csv"
    for p in (mcw, atlas_i, atlas_ii):
        if not p.exists():
            raise FileNotFoundError(f"missing {p} — run: just download")

    frames: list[pl.DataFrame] = []
    protein_only_rows: list[pl.DataFrame] = []

    odb = OGlcnacDbAdapter()
    obs_db, po_db = odb.parse(mcw)
    frames.append(obs_db)
    protein_only_rows.append(po_db)
    log.info("oglcnacdb_parsed", rows=obs_db.height, protein_only=po_db.height)

    for adapter_cls, path in [
        (AtlasDatasetIAdapter, atlas_i),
        (AtlasDatasetIIAdapter, atlas_ii),
    ]:
        adapter = adapter_cls()
        obs = adapter.parse(path)
        frames.append(obs)
        log.info("atlas_parsed", dataset=adapter.dataset_id, rows=obs.height)

    obs_all = pl.concat(frames, how="diagonal_relaxed")
    obs_staging = MT / "staging"
    obs_staging.mkdir(parents=True, exist_ok=True)
    obs_all.write_parquet(obs_staging / "canon_observations.parquet")

    resolved, identity_dim = resolve_identity(obs_all)

    protein_only_count = sum(po.height for po in protein_only_rows)

    encoded, dims = encode_context(resolved)
    with_uid = assign_set_uid(encoded)
    write_registry(identity_dim, dims)

    for dataset_id in with_uid["dataset_id"].unique().to_list():
        part = with_uid.filter(pl.col("dataset_id") == dataset_id)
        coords = part
        metrics = part
        write_partition(coords, metrics, dataset_id)

    FIGURES.mkdir(parents=True, exist_ok=True)
    completeness = axis_completeness(with_uid)
    completeness.write_csv(FIGURES / "axis_completeness_canon.csv")
    overlap = canon_overlap(with_uid)
    overlap.write_json(FIGURES / "canon_overlap.json")

    summary = {
        "observation_rows": obs_all.height,
        "unique_sites": identity_dim.filter(~pl.col("protein_level_only")).height,
        "unique_sets": with_uid["set_uid"].n_unique(),
        "datasets": with_uid["dataset_id"].unique().to_list(),
        "protein_level_only": protein_only_count,
    }
    (MT / "canon_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("canon_complete", **summary)
    print(json.dumps(summary, indent=2))
