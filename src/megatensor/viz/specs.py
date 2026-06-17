"""Validated figure wrappers — agents fill specs, not matplotlib knobs (FIGURES.md §6)."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from megatensor.viz import preamble  # noqa: F401 — side-effect: Agg + rcParams
from megatensor.viz.theme import apply_seaborn_theme, despine

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

apply_seaborn_theme()

_PANEL_W = 11.0
_PANEL_H_SM = 4.5
_PANEL_H_MD = 6.0
_PANEL_H_LG = 7.5


def _apply_panel_theme() -> None:
    apply_seaborn_theme()


def _axes_list(axes) -> list:
    if isinstance(axes, np.ndarray):
        return list(axes.flat)
    if isinstance(axes, dict):
        return [axes[k] for k in sorted(axes)]
    return [axes]


def _finish_bar_ax(ax) -> None:
    despine(ax)
    ax.yaxis.grid(True, color="#ebebeb", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def _finish_cat_ax(ax, *, horizontal: bool = False) -> None:
    _finish_bar_ax(ax)
    if horizontal:
        ax.tick_params(axis="y", length=0)
    else:
        ax.tick_params(axis="x", length=0)


def _scale_label(values: Sequence[float], base: str) -> tuple[np.ndarray, str]:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or np.nanmax(np.abs(arr)) < 1e4:
        return arr, base
    if np.nanmax(np.abs(arr)) >= 1e6:
        return arr / 1e6, f"{base} (millions)"
    if np.nanmax(np.abs(arr)) >= 1e3:
        return arr / 1e3, f"{base} (thousands)"
    return arr, base


def save_figure(fig: Figure, stem: Path) -> dict[str, str]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}
    pdf = stem.with_suffix(".pdf")
    fig.savefig(pdf, metadata={"Date": None}, facecolor="white", bbox_inches="tight")
    paths["pdf"] = str(pdf)
    png = stem.with_suffix(".png")
    fig.savefig(png, facecolor="white", bbox_inches="tight")
    paths["png"] = str(png)
    return paths


def completeness_axis_profile(
    canon_matrix: np.ndarray,
    pride_matrix: np.ndarray,
    col_labels: Sequence[str],
    *,
    title: str,
) -> Figure:
    """Figure A: mean axis fill % — canon vs PRIDE (no heatmap)."""
    canon_mean = np.nanmean(canon_matrix, axis=0)
    pride_mean = np.nanmean(pride_matrix, axis=0)
    rows = []
    for lab, c, p in zip(col_labels, canon_mean, pride_mean, strict=True):
        rows.append({"axis": lab, "pct": c, "layer": "Canon"})
        rows.append({"axis": lab, "pct": p, "layer": "PRIDE"})
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(_PANEL_W, _PANEL_H_SM), constrained_layout=True)
    sns.barplot(
        data=df,
        x="axis",
        y="pct",
        hue="layer",
        palette={"Canon": preamble.COLOR_CANON, "PRIDE": preamble.COLOR_PRIDE},
        ax=ax,
        edgecolor="white",
        width=0.72,
    )
    ax.set_ylabel("Mean rows with axis filled (%)")
    ax.set_xlabel("")
    ax.set_ylim(0, 105)
    ax.legend(title="Layer", frameon=False, loc="upper right")
    ax.set_title(title, pad=12)
    _finish_cat_ax(ax)
    return fig


def completeness_heatmap_split(
    canon_matrix: np.ndarray,
    canon_rows: Sequence[str],
    pride_matrix: np.ndarray,
    pride_rows: Sequence[str],
    col_labels: Sequence[str],
    *,
    title: str,
) -> Figure:
    """Deprecated alias — bar profile only."""
    return completeness_axis_profile(canon_matrix, pride_matrix, col_labels, title=title)


def overlap_bars(
    panels: Sequence[tuple[dict[str, int], str]],
    *,
    title: str,
) -> Figure:
    """Figure B: shared sites highlighted; % of union annotated."""
    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(_PANEL_W, _PANEL_H_SM), sharey=True, constrained_layout=True)
    axes_list = _axes_list(axes) if n > 1 else [axes]

    keys = ["canon_only", "shared_sites", "pride_only"]
    tick_labels = ["Canon only", "Shared", "PRIDE only"]
    ymax = 0

    for ax, (stats, subtitle) in zip(axes_list, panels, strict=True):
        vals = [stats[k] for k in keys]
        ymax = max(ymax, max(vals))
        total = sum(vals)
        colors = [preamble.COLOR_MUTED, preamble.COLOR_HIGHLIGHT, preamble.COLOR_MUTED]
        df = pd.DataFrame({"category": tick_labels, "sites": vals, "color": colors})
        bars = sns.barplot(
            data=df,
            x="category",
            y="sites",
            hue="category",
            palette=dict(zip(tick_labels, colors, strict=True)),
            dodge=False,
            legend=False,
            ax=ax,
            width=0.62,
            edgecolor="white",
            linewidth=1.2,
            zorder=3,
        )
        ax.set_title(subtitle, fontsize=13, pad=10)
        _finish_cat_ax(ax)
        ax.set_xlim(-0.6, len(tick_labels) - 0.4)

    axes_list[0].set_ylabel("Unique sites")
    for ax in axes_list[1:]:
        ax.tick_params(axis="y", labelleft=False)
    for ax in axes_list:
        ax.set_ylim(0, ymax * 1.18)
    fig.suptitle(title, fontsize=16, fontweight="bold", y=1.05)
    return fig


def ranked_hbar(
    labels: Sequence[str],
    values: Sequence[float],
    *,
    title: str,
    xlabel: str,
    engines: Sequence[str] | None = None,
) -> Figure:
    """Horizontal ranked bars — engine encoded by color, legend outside plot."""
    order = np.argsort(values)
    labels = [labels[i] for i in order]
    values = [values[i] for i in order]
    engine_labels = None
    if engines is not None:
        engines = [engines[i] for i in order]
        engine_labels = [
            "MaxQuant" if "maxquant" in str(e).lower() or "mq" in str(e).lower() else "Proteome Discoverer"
            for e in engines
        ]

    height = max(_PANEL_H_MD, 0.42 * len(labels) + 1.0)
    fig, ax = plt.subplots(figsize=(9.5, height), constrained_layout=True)
    df = pd.DataFrame({"pxd": labels, "sites": values})
    if engine_labels:
        df["engine"] = engine_labels
        sns.barplot(
            data=df,
            y="pxd",
            x="sites",
            hue="engine",
            dodge=False,
            palette={"MaxQuant": preamble.COLOR_MQ, "Proteome Discoverer": preamble.COLOR_PD},
            ax=ax,
            edgecolor="white",
            zorder=3,
        )
        ax.legend(title="Engine", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    else:
        sns.barplot(data=df, y="pxd", x="sites", color=preamble.COLOR_MUTED, ax=ax, edgecolor="white", zorder=3)

    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("")
    ax.set_title(title, pad=12)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    _finish_cat_ax(ax, horizontal=True)
    return fig


def condition_bar(
    categories: Sequence[str],
    values: Sequence[float],
    counts: Sequence[int],
    *,
    title: str,
    subtitle: str = "",
    ylabel_base: str = "Mean intensity",
) -> Figure:
    """Figure C: SILAC comparison — metadata in legend, not on bars."""
    scaled, ylab = _scale_label(values, ylabel_base)
    color_map = {"light": "#56B4E9", "heavy": "#D55E00"}
    legend_labels = []
    for cat, n in zip(categories, counts, strict=True):
        legend_labels.append(f"{str(cat).capitalize()} (n={n})")

    fc_note = ""
    if len(scaled) == 2 and scaled[0] > 0 and scaled[1] > 0:
        heavy_i = 0 if str(categories[0]).lower() == "heavy" else 1
        light_i = 1 - heavy_i
        fc_note = f"{scaled[heavy_i] / scaled[light_i]:.2f}× Heavy/Light"

    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    df = pd.DataFrame(
        {
            "category": list(categories),
            "intensity": scaled,
            "hue": legend_labels,
        }
    )
    palette = {lab: color_map.get(str(cat).lower(), preamble.COLOR_MUTED) for cat, lab in zip(categories, legend_labels, strict=True)}
    sns.barplot(
        data=df,
        x="category",
        y="intensity",
        hue="hue",
        palette=palette,
        dodge=False,
        ax=ax,
        width=0.55,
        edgecolor="white",
        linewidth=1.2,
        zorder=3,
    )
    ax.set_ylabel(ylab)
    ax.set_xlabel("")
    ax.set_title(title, pad=14)
    ax.set_ylim(0, float(np.nanmax(scaled)) * 1.12 if len(scaled) else 1)

    legend_title = " · ".join(p for p in (fc_note, subtitle) if p)
    ax.legend(title=legend_title or None, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    _finish_cat_ax(ax)
    return fig


def feature_coverage_hbar(
    labels: Sequence[str],
    percents: Sequence[float],
    *,
    title: str,
) -> Figure:
    """Enrichment completeness — drop zero rows, cap axis to data."""
    pairs = [(lab, p) for lab, p in zip(labels, percents, strict=True) if p > 0.5]
    if not pairs:
        pairs = list(zip(labels, percents, strict=True))
    labels, percents = zip(*pairs, strict=True) if pairs else ([], [])

    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.9 * len(labels) + 1.5)), constrained_layout=True)
    colors = [preamble.COLOR_HIGHLIGHT if p >= 80 else "#7eb6d9" if p >= 60 else preamble.COLOR_MUTED for p in percents]
    df = pd.DataFrame({"label": labels, "pct": percents, "color": colors})
    bars = sns.barplot(
        data=df,
        y="label",
        x="pct",
        hue="label",
        palette=dict(zip(labels, colors, strict=True)),
        dodge=False,
        legend=False,
        ax=ax,
        edgecolor="white",
        linewidth=1.0,
        zorder=3,
    )
    ax.invert_yaxis()
    xmax = min(100, max(percents) * 1.12) if len(percents) else 100
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Sites with feature (%)")
    ax.set_ylabel("")
    ax.set_title(title, pad=12)
    _finish_cat_ax(ax, horizontal=True)
    return fig
