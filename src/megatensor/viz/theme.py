"""Seaborn theme + adjustText label repulsion (FIGURES.md)."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from megatensor.viz import preamble

import matplotlib.pyplot as plt
import seaborn as sns


def apply_seaborn_theme() -> None:
    """House style via mplstyle + seaborn whitegrid panel theme."""
    plt.style.use(preamble.STYLE_PATH)
    sns.set_theme(
        style="whitegrid",
        context="notebook",
        font_scale=1.05,
        rc={
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.0,
            "axes.titlesize": 15,
            "axes.titleweight": "bold",
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 10,
            "legend.fontsize": 11,
            "figure.dpi": 150,
            "savefig.dpi": 200,
            "grid.color": "#ebebeb",
            "grid.linewidth": 0.7,
            "axes.axisbelow": True,
        },
    )


def despine(ax, *, left: bool = True, bottom: bool = True) -> None:
    sns.despine(ax=ax, top=True, right=True, left=not left, bottom=not bottom, trim=False)


def _cmap(name: str):
    try:
        import cmcrameri.cm as cmc

        return getattr(cmc, name)
    except Exception:
        return sns.color_palette(name, as_cmap=True)


def matrix_heatmap(
    data: np.ndarray,
    *,
    ax=None,
    xticklabels=True,
    yticklabels=True,
    cmap: str = "batlow",
    vmin=None,
    vmax=None,
    center=None,
    mask: np.ndarray | None = None,
    annot: bool = False,
    fmt: str = ".2f",
    square: bool = False,
    cbar_label: str = "",
    yticklabelsize: float = 8,
    xticklabelsize: float = 9,
):
    """Clean heatmap — no grid lines, perceptually uniform cmap."""
    cmap_obj = _cmap(cmap) if isinstance(cmap, str) else cmap
    hm = sns.heatmap(
        data,
        ax=ax,
        mask=mask,
        cmap=cmap_obj,
        vmin=vmin,
        vmax=vmax,
        center=center,
        square=square,
        annot=annot,
        fmt=fmt,
        annot_kws={"size": 9, "color": "0.15"},
        linewidths=0,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        cbar_kws={"label": cbar_label, "shrink": 0.82, "aspect": 28},
        rasterized=True,
    )
    if ax is not None and yticklabels is not False:
        ax.tick_params(axis="y", labelsize=yticklabelsize)
        ax.tick_params(axis="x", labelsize=xticklabelsize)
    return hm


def cluster_heatmap(
    data: pd.DataFrame,
    *,
    cmap: str = "batlow",
    vmin=None,
    vmax=None,
    cbar_label: str = "",
    figsize: tuple[float, float] = (10, 8),
    row_cluster: bool = True,
    col_cluster: bool = False,
) -> sns.matrix.ClusterGrid:
    """Clustered heatmap for site × feature matrices."""
    plot_df = data.astype(float).fillna(0.0)
    cmap_obj = _cmap(cmap)
    with plt.rc_context({"figure.constrained_layout.use": False}):
        g = sns.clustermap(
            plot_df,
            cmap=cmap_obj,
            vmin=vmin,
            vmax=vmax,
            row_cluster=row_cluster and plot_df.shape[0] > 1,
            col_cluster=col_cluster and plot_df.shape[1] > 1,
            linewidths=0,
            figsize=figsize,
            dendrogram_ratio=0.08,
            cbar_kws={"label": cbar_label},
            yticklabels=True,
            xticklabels=True,
            rasterized=True,
        )
    g.fig.set_constrained_layout(False)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=7)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=9, rotation=45, ha="right")
    g.fig.subplots_adjust(top=0.92)
    return g


def repel_labels(
    ax,
    x: Sequence[float],
    y: Sequence[float],
    labels: Sequence[str],
    *,
    fontsize: float = 8,
    color: str = "0.15",
    scatter_x: Sequence[float] | None = None,
    scatter_y: Sequence[float] | None = None,
    ha: str = "left",
    va: str = "center",
    force_text: tuple[float, float] = (0.2, 0.4),
    force_points: tuple[float, float] = (0.1, 0.2),
    expand_points: tuple[float, float] = (1.2, 1.4),
    arrowprops: dict | None = None,
) -> list:
    """Place labels with adjustText; optional scatter repulsion field."""
    from adjustText import adjust_text

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    texts = [
        ax.text(float(xi), float(yi), str(lab), fontsize=fontsize, color=color, ha=ha, va=va)
        for xi, yi, lab in zip(x, y, labels, strict=True)
    ]
    sx = np.asarray(scatter_x if scatter_x is not None else x, dtype=float)
    sy = np.asarray(scatter_y if scatter_y is not None else y, dtype=float)
    adjust_text(
        texts,
        x=sx,
        y=sy,
        ax=ax,
        force_text=force_text,
        force_points=force_points,
        expand_points=expand_points,
        arrowprops=arrowprops or {"arrowstyle": "-", "color": "0.55", "lw": 0.6, "alpha": 0.8, "shrinkA": 4, "shrinkB": 2},
    )
    return texts


def bar_value_labels(
    ax,
    bars,
    labels: Sequence[str],
    *,
    fontsize: float = 9,
    padding: float = 0.01,
    horizontal: bool = False,
    repel: bool = True,
) -> None:
    """Annotate bar ends; use adjustText when repel=True."""
    xs: list[float] = []
    ys: list[float] = []
    for bar, lab in zip(bars, labels, strict=True):
        if horizontal:
            xs.append(bar.get_width() + padding * (ax.get_xlim()[1] - ax.get_xlim()[0]))
            ys.append(bar.get_y() + bar.get_height() / 2)
        else:
            xs.append(bar.get_x() + bar.get_width() / 2)
            ys.append(bar.get_height() + padding * (ax.get_ylim()[1] - ax.get_ylim()[0]))

    if repel and len(labels) > 1:
        repel_labels(ax, xs, ys, labels, fontsize=fontsize, ha="left" if horizontal else "center", va="center")
    else:
        for xi, yi, lab in zip(xs, ys, labels, strict=True):
            ax.text(xi, yi, lab, fontsize=fontsize, ha="left" if horizontal else "center", va="bottom" if not horizontal else "center")
