"""Figure A/B data products and trajectory candidates."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.paths import FIGURES
from megatensor.pipeline.union import run_union
from megatensor.store import CANON_STORE, PRIDE_STORE
from megatensor.viz.upset import export_upset_tables

try:
    from megatensor.viz.plots import render_panel_figures
except ImportError:
    render_panel_figures = None  # type: ignore[misc, assignment]

log = structlog.get_logger()


def run_figures() -> dict:
    """Emit separated figure datasets for canon vs PRIDE narratives."""
    if not CANON_STORE.summary_path.exists():
        raise FileNotFoundError("run: just canon")
    if not PRIDE_STORE.summary_path.exists():
        raise FileNotFoundError("run: just pride-tensorize")

    FIGURES.mkdir(parents=True, exist_ok=True)
    union_summary = run_union()

    # Figure C candidate: conditioned PRIDE site with multi-metric spread (Light/Heavy)
    traj = _trajectory_candidate()
    if traj is not None:
        traj.write_csv(FIGURES / "figure_c_trajectory_candidate.csv")

    upset = export_upset_tables()
    upset["canon_overlap_pairs"].write_csv(FIGURES / "upset_canon_pairwise.csv")
    upset["canon_membership"].write_csv(FIGURES / "upset_canon_membership.csv")
    upset["pride_membership"].write_csv(FIGURES / "upset_pride_membership.csv")

    out = {
        "canon_summary": json.loads(CANON_STORE.summary_path.read_text()),
        "pride_summary": json.loads(PRIDE_STORE.summary_path.read_text()),
        "union_summary": union_summary,
        "figure_c_rows": traj.height if traj is not None else 0,
    }
    if render_panel_figures is not None:
        try:
            out["panel_figures"] = render_panel_figures()
        except Exception as exc:
            log.warning("panel_figures_failed", error=str(exc))
    (FIGURES / "figures_manifest.json").write_text(json.dumps(out, indent=2))
    log.info("figures_complete", **{k: v for k, v in out.items() if k != "canon_summary"})
    print(json.dumps(out, indent=2))
    return out


def _trajectory_candidate() -> pl.DataFrame | None:
    """PXD039536 Light/Heavy site tables — multi-condition on same sites."""
    obs_path = PRIDE_STORE.staging / "observations.parquet"
    if not obs_path.exists():
        return None
    obs = pl.read_parquet(obs_path)
    cand = obs.filter(pl.col("dataset_id") == "PXD039536")
    if cand.is_empty():
        return None
    top = (
        cand.group_by(["protein_id_raw", "residue_pos_raw", "residue_aa"])
        .agg(
            pl.col("cond_treatment").drop_nulls().unique().alias("treatments"),
            pl.col("metric_value").mean().alias("mean_intensity"),
            pl.len().alias("n_obs"),
        )
        .filter(pl.col("treatments").list.len() >= 2)
        .sort("n_obs", descending=True)
        .head(20)
    )
    if top.is_empty():
        return cand.head(50)
    keys = top.head(1)
    site = keys.row(0, named=True)
    return cand.filter(
        (pl.col("protein_id_raw") == site["protein_id_raw"])
        & (pl.col("residue_pos_raw") == site["residue_pos_raw"])
        & (pl.col("residue_aa") == site["residue_aa"])
    ).select(
        "dataset_id",
        "protein_id_raw",
        "residue_pos_raw",
        "residue_aa",
        "cond_treatment",
        "metric_name",
        "metric_value",
        "loc_score",
        "inst_model",
        "prov_country",
    )
