"""Phase 0: canon reference libraries -> isolated canon tensor."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.ingest.canon_atlas import AtlasDatasetIAdapter, AtlasDatasetIIAdapter
from megatensor.ingest.canon_oglcnacdb import OGlcnacDbAdapter
from megatensor.paths import CANON, FIGURES
from megatensor.pipeline.tensorize import tensorize_observations
from megatensor.store import CANON_STORE
from megatensor.viz.completeness import canon_overlap

log = structlog.get_logger()


def run_canon() -> dict:
    CANON.mkdir(parents=True, exist_ok=True)

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
    protein_only_count = sum(po.height for po in protein_only_rows)

    summary, with_uid = tensorize_observations(
        obs_all,
        CANON_STORE,
        extra_summary={
            "purpose": "reference_library_harmonization",
            "protein_level_only": protein_only_count,
        },
    )

    FIGURES.mkdir(parents=True, exist_ok=True)
    canon_overlap(with_uid).write_json(FIGURES / "canon_vs_canon_overlap.json")

    print(json.dumps(summary, indent=2))
    return summary
