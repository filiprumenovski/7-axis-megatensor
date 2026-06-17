"""Panel figure entrypoints — thin adapters over viz/specs wrappers."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl

from megatensor.paths import FIGURES
from megatensor.viz import preamble
from megatensor.viz.specs import (
    completeness_heatmap_split,
    condition_bar,
    feature_coverage_hbar,
    overlap_bars,
    ranked_hbar,
    save_figure,
)

AXIS_ORDER = [
    "identity",
    "ptm",
    "quant",
    "condition",
    "acquisition",
    "instrument",
    "provenance",
]
AXIS_LABELS = {
    "identity": "Identity",
    "ptm": "PTM",
    "quant": "Quant",
    "condition": "Condition",
    "acquisition": "Acquisition",
    "instrument": "Instrument",
    "provenance": "Provenance",
}


def _export(fig, stem: str) -> dict[str, str]:
    paths = save_figure(fig, FIGURES / stem)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return paths


def plot_axis_completeness() -> dict[str, str]:
    canon_path = FIGURES / "axis_completeness_canon.csv"
    pride_path = FIGURES / "axis_completeness_pride.csv"
    if not canon_path.is_file() or not pride_path.is_file():
        raise FileNotFoundError("run tensorize first — axis completeness CSVs missing")

    canon = pl.read_csv(canon_path)
    pride = pl.read_csv(pride_path)
    col_labels = [AXIS_LABELS[a] for a in AXIS_ORDER]
    cols = [f"pct_{a}" for a in AXIS_ORDER]

    fig = completeness_heatmap_split(
        canon.select(cols).to_numpy(),
        canon["dataset_id"].to_list(),
        pride.select(cols).to_numpy(),
        pride["dataset_id"].to_list(),
        col_labels,
        title="Canon fills identity; PRIDE fills instrument and provenance",
    )
    return _export(fig, "figure_a_axis_completeness")


def plot_cross_layer_overlap() -> dict[str, str]:
    overlap_path = FIGURES / "canon_vs_pride_overlap.json"
    if not overlap_path.is_file():
        raise FileNotFoundError("run union first")

    data = {row["metric"]: int(row["value"]) for row in json.loads(overlap_path.read_text())}
    panels: list[tuple[dict[str, int], str]] = [(data, "All PRIDE deposits")]
    human_path = FIGURES / "canon_vs_pride_overlap_human.json"
    if human_path.is_file():
        human = {row["metric"]: int(row["value"]) for row in json.loads(human_path.read_text())}
        panels.append((human, "Human PRIDE (rice excluded)"))

    fig = overlap_bars(
        panels,
        title=f"{data['shared_sites']:,} sites appear in both canon and PRIDE",
    )
    return _export(fig, "figure_b_cross_layer_overlap")


def plot_pride_heterogeneity() -> dict[str, str]:
    spread_path = FIGURES / "pride_engine_spread.csv"
    if not spread_path.is_file():
        raise FileNotFoundError("run pride-tensorize first")

    df = pl.read_csv(spread_path).sort("unique_sites", descending=False)
    fig = ranked_hbar(
        df["dataset_id"].to_list(),
        df["unique_sites"].to_list(),
        title="PRIDE deposit site counts",
        xlabel="Unique O-GlcNAc sites",
        engines=df["engines"].to_list(),
    )
    return _export(fig, "figure_pride_heterogeneity")


def plot_trajectory() -> dict[str, str]:
    traj_path = FIGURES / "figure_c_trajectory_candidate.csv"
    if not traj_path.is_file():
        raise FileNotFoundError("run figures first")

    df = pl.read_csv(traj_path)
    if df.is_empty():
        raise RuntimeError("trajectory candidate empty")

    site = df.row(0, named=True)
    site_label = f"{site['protein_id_raw']}:{site['residue_pos_raw']}{site['residue_aa']}"
    pxd = site["dataset_id"]

    inten = df.filter(pl.col("metric_name") == "intensity")
    if inten.is_empty():
        inten = df
    summary = (
        inten.group_by("cond_treatment")
        .agg(pl.col("metric_value").mean().alias("mean_intensity"), pl.len().alias("n"))
        .sort("cond_treatment")
    )

    fig = condition_bar(
        summary["cond_treatment"].to_list(),
        summary["mean_intensity"].to_list(),
        summary["n"].to_list(),
        title=f"{site_label} is heavier under SILAC",
        subtitle=f"{pxd} · Orbitrap Fusion Lumos · China",
        ylabel_base="Mean intensity",
    )
    return _export(fig, "figure_c_trajectory")


def plot_enrichment_completeness() -> dict[str, str]:
    comp_path = FIGURES / "enrichment_completeness.json"
    if not comp_path.is_file():
        raise FileNotFoundError("run enrich first")

    data = json.loads(comp_path.read_text())
    keys = ["gene_symbol_pct", "domain_pct", "seq_window_pct", "disorder_pct"]
    labels = ["Gene symbol", "Domain / region", "Sequence window", "Disorder (metapredict)"]
    vals = [float(data.get(k, 0)) for k in keys]

    fig = feature_coverage_hbar(
        labels,
        vals,
        title=f"UniProt features cover ~{np.nanmax(vals):.0f}% of union sites",
    )
    return _export(fig, "figure_enrichment_completeness")


def render_panel_figures() -> dict[str, dict[str, str]]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    builders = [
        ("figure_a", plot_axis_completeness),
        ("figure_b", plot_cross_layer_overlap),
        ("figure_c", plot_trajectory),
        ("pride_heterogeneity", plot_pride_heterogeneity),
    ]
    paths: dict[str, dict[str, str]] = {}
    for key, fn in builders:
        paths[key] = fn()

    comp_path = FIGURES / "enrichment_completeness.json"
    if comp_path.is_file():
        paths["enrichment"] = plot_enrichment_completeness()

    (FIGURES / "panel_figures.json").write_text(json.dumps(paths, indent=2))
    return paths
