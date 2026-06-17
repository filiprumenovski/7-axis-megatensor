"""Route downloaded PRIDE files to engine-specific adapters."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from megatensor.ingest.pride_collapse import collapse_pride_observations
from megatensor.ingest.pride_common import extract_from_zip, finalize_obs
from megatensor.ingest.pride_diann import parse_diann_zip
from megatensor.ingest.pride_mq_sites import parse_mq_sites
from megatensor.ingest.pride_mztab import parse_mztab
from megatensor.ingest.pride_pd import parse_pd_psm
from megatensor.ingest.pride_xlsx import parse_pd_xlsx
from megatensor.ingest.pride_download import build_manifest
from megatensor.paths import PRIDE_DOWNLOADS


def parse_file(path: Path, engine_kind: str, pxd: str) -> pl.DataFrame:
    if engine_kind in ("maxquant_sites", "site_table"):
        return parse_mq_sites(path, pxd=pxd)
    if engine_kind == "site_zip":
        inner = extract_from_zip(path)
        return parse_mq_sites(inner, pxd=pxd)
    if engine_kind == "pd_psm":
        return parse_pd_psm(path, pxd=pxd)
    if engine_kind == "mztab_result":
        return parse_mztab(path, pxd=pxd)
    if engine_kind == "diann_zip":
        return parse_diann_zip(path, pxd=pxd)
    if engine_kind == "skyline_reporter":
        # Scan-level reporter ions without peptide/site mapping — not SET-able.
        return finalize_obs(pl.DataFrame())
    if engine_kind == "csv_table":
        return finalize_obs(pl.DataFrame())
    if engine_kind == "xlsx_table":
        return parse_pd_xlsx(path, pxd=pxd)
    return finalize_obs(pl.DataFrame())


def parse_downloads(download_root: Path | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    root = download_root or PRIDE_DOWNLOADS
    manifest = build_manifest()
    frames: list[pl.DataFrame] = []
    qc_rows: list[dict] = []

    for row in manifest.iter_rows(named=True):
        pxd = row["project_accession"]
        fname = row["file_name"]
        path = root / pxd / fname
        engine = row["engine_kind"]
        if not path.is_file():
            qc_rows.append(
                {
                    "pxd": pxd,
                    "file_name": fname,
                    "engine_kind": engine,
                    "status": "missing",
                    "observation_rows": 0,
                    "unique_sites": 0,
                }
            )
            continue
        try:
            obs = parse_file(path, engine, pxd)
        except Exception as exc:
            qc_rows.append(
                {
                    "pxd": pxd,
                    "file_name": fname,
                    "engine_kind": engine,
                    "status": f"error:{exc}",
                    "observation_rows": 0,
                    "unique_sites": 0,
                }
            )
            continue
        sites = 0
        if obs.height:
            sites = obs.select(["protein_id_raw", "residue_pos_raw", "residue_aa"]).unique().height
        qc_rows.append(
            {
                "pxd": pxd,
                "file_name": fname,
                "engine_kind": engine,
                "status": "ok" if obs.height else "empty",
                "observation_rows": obs.height,
                "unique_sites": sites,
            }
        )
        if obs.height:
            frames.append(obs)

    obs_all = pl.concat(frames, how="diagonal_relaxed") if frames else finalize_obs(pl.DataFrame())
    obs_all = collapse_pride_observations(obs_all)
    qc = pl.DataFrame(qc_rows)
    return obs_all, qc
