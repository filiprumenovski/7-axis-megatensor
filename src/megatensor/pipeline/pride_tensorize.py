"""Phase 3: PRIDE experimental deposits -> isolated pride tensor."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.ingest.pride_router import parse_downloads
from megatensor.paths import FIGURES, PRIDE_DOWNLOADS
from megatensor.pipeline.tensorize import tensorize_observations
from megatensor.store import PRIDE_STORE
from megatensor.viz.completeness import pride_engine_spread

log = structlog.get_logger()


def run_pride_tensorize(download_root=None, *, reparse: bool = False) -> dict:
    from pathlib import Path

    root = Path(download_root) if download_root else PRIDE_DOWNLOADS
    obs_path = PRIDE_STORE.staging / "observations.parquet"

    if reparse or not obs_path.exists():
        if not root.exists():
            raise FileNotFoundError(f"No downloads at {root} — run: just pride-download")
        obs, qc = parse_downloads(root)
        PRIDE_STORE.staging.mkdir(parents=True, exist_ok=True)
        obs.write_parquet(obs_path)
        qc.write_csv(PRIDE_STORE.staging / "parse_qc.csv")
        FIGURES.mkdir(parents=True, exist_ok=True)
        qc.write_csv(FIGURES / "pride_ingest_files.csv")
    else:
        obs = pl.read_parquet(obs_path)

    if obs.is_empty():
        raise RuntimeError("PRIDE observations empty — check adapters / downloads")

    summary, with_uid = tensorize_observations(
        obs,
        PRIDE_STORE,
        extra_summary={
            "purpose": "experimental_deposit_tensorization",
            "pxds": obs["dataset_id"].unique().sort().to_list(),
        },
    )

    pride_engine_spread(with_uid).write_csv(FIGURES / "pride_engine_spread.csv")
    print(json.dumps(summary, indent=2))
    return summary
