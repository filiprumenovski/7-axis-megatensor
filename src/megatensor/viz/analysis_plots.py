"""Figures driven by megatensor/analysis outputs — findings, not QC."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from megatensor.paths import FIGURES, ROOT
from megatensor.viz import preamble
from megatensor.viz.specs import _finish_bar_ax, _finish_cat_ax, save_figure
from megatensor.viz.theme import apply_seaborn_theme, repel_labels

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

apply_seaborn_theme()

ANALYSIS = ROOT / "megatensor" / "analysis"
_PANEL_W = 11.0


def _export(fig, stem: str) -> dict[str, str]:
    paths = save_figure(fig, FIGURES / stem)
    plt.close(fig)
    return paths


def _load(name: str) -> pl.DataFrame:
    p = ANALYSIS / f"{name}.parquet"
    if not p.is_file():
        raise FileNotFoundError(f"missing {p} — run: megatensor analyze")
    return pl.read_parquet(p)



def plot_silac_ma(*, label_top: int = 8) -> dict[str, str]:
    """SILAC M–A: log2(Heavy/Light) vs mean abundance; label top-|M| sites."""
    fc = _load("silac_fc")
    df = fc.to_pandas()
    m = df["log2_fc"].to_numpy()
    a = df["log2_a"].to_numpy()
    names = df["site_key"].tolist()

    fig, ax = plt.subplots(figsize=(_PANEL_W, 6), constrained_layout=True)
    sns.scatterplot(data=df, x="log2_a", y="log2_fc", ax=ax, s=18, color=preamble.COLOR_MUTED, alpha=0.55, edgecolor=None, rasterized=True)

    order = np.argsort(np.abs(m))[::-1][:label_top]
    hi = df.iloc[order]
    sns.scatterplot(
        data=hi,
        x="log2_a",
        y="log2_fc",
        ax=ax,
        s=48,
        color=preamble.COLOR_HIGHLIGHT,
        edgecolor="white",
        linewidth=0.5,
        zorder=4,
    )
    repel_labels(
        ax,
        hi["log2_a"].to_numpy(),
        hi["log2_fc"].to_numpy(),
        hi["site_key"].tolist(),
        scatter_x=a,
        scatter_y=m,
        fontsize=8,
    )

    med = float(np.median(m))
    ax.axhline(0, color="#888888", lw=0.8, ls="--")
    ax.axhline(med, color=preamble.COLOR_ACCENT, lw=1.0, ls=":", label=f"median M={med:.2f}")
    ax.set_xlabel("log₂ A (mean abundance)")
    ax.set_ylabel("log₂(Heavy / Light)")
    ax.set_title(f"PXD039536 SILAC — {fc.height} sites; median {np.median(np.exp2(m)):.2f}× Heavy/Light")
    _finish_bar_ax(ax)
    ax.legend(loc="upper right", frameon=False)
    return _export(fig, "analysis_silac_ma")


def plot_concordance_scatter() -> dict[str, str]:
    pair = _load("concordance_silac_pair")
    summary = json.loads((FIGURES / "analysis_summary.json").read_text())
    stats = summary.get("concordance", {})
    pdf = pair.to_pandas()

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    sns.scatterplot(data=pdf, x="log10_a", y="log10_b", ax=ax, s=28, color=preamble.COLOR_HIGHLIGHT, alpha=0.65, edgecolor="white", linewidth=0.4)
    x, y = pdf["log10_a"].to_numpy(), pdf["log10_b"].to_numpy()
    if len(x) > 2:
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, coef[0] * xs + coef[1], color=preamble.COLOR_NEG, lw=1.2, ls="--", label="OLS fit")
    r = stats.get("pearson_r_log10", 0)
    n = stats.get("n_shared_sites", pair.height)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="w", markerfacecolor=preamble.COLOR_HIGHLIGHT, markersize=8, label=f"Shared sites (n={n})"),
            Line2D([0], [0], color=preamble.COLOR_NEG, ls="--", label=f"Pearson r = {r}"),
        ],
        frameon=False,
        loc="lower right",
    )
    ax.set_xlabel("PXD039536 (China SILAC)")
    ax.set_ylabel("PXD058744 (US MaxQuant)")
    ax.set_title("Shared-site intensities correlate weakly across studies")
    _finish_bar_ax(ax)
    return _export(fig, "analysis_concordance_scatter")


def plot_replication_tiers() -> dict[str, str]:
    tiers = _load("overlap_tiers")
    pretty = {
        "pride_unique": "PRIDE-only",
        "pride_multi_pxd": "PRIDE ≥2 PXDs",
        "canon_and_pride": "Canon ∩ PRIDE",
        "canon_and_pride_multi_pxd": "Triangulated\n(canon + ≥2 PXDs)",
    }
    colors = [preamble.COLOR_MUTED, "#7eb6d9", "#7eb6d9", preamble.COLOR_HIGHLIGHT]
    df = tiers.to_pandas()
    df["label"] = df["tier"].map(pretty)
    df["color"] = colors

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    bars = sns.barplot(
        data=df,
        x="label",
        y="n_sites",
        hue="label",
        palette=dict(zip(df["label"], colors, strict=True)),
        dodge=False,
        legend=False,
        ax=ax,
        width=0.62,
        edgecolor="white",
        linewidth=1.2,
    )
    ax.set_ylabel("Unique sites")
    ax.set_xlabel("")
    ax.set_title("Replication tiers in the O-GlcNAc megatensor")
    _finish_cat_ax(ax)
    return _export(fig, "analysis_replication_tiers")


def plot_protein_hubs(*, top_n: int = 12) -> dict[str, str]:
    hubs = _load("protein_hubs").head(top_n)
    rows = []
    for r in hubs.iter_rows(named=True):
        gene = r.get("gene_symbol") or "?"
        rows.append(
            {
                "label": f"{gene} ({r['protein_acc']}) · {r['max_pxds']} PXDs",
                "n": r["n_triangulated_sites"],
            }
        )
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(_PANEL_W, 5), constrained_layout=True)
    sns.barplot(
        data=df,
        y="label",
        x="n",
        color=preamble.COLOR_HIGHLIGHT,
        edgecolor="white",
        ax=ax,
        zorder=3,
    )
    ax.invert_yaxis()
    ax.set_xlabel("Triangulated sites (canon + multi-PXD)")
    ax.set_ylabel("")
    ax.set_title("Proteins with the most cross-study O-GlcNAc evidence")
    _finish_cat_ax(ax, horizontal=True)
    return _export(fig, "analysis_protein_hubs")


def plot_chemoproteomics(*, top_n: int = 12) -> dict[str, str]:
    mat = _load("chemoproteomics_matrix")
    probe_cols = [c for c in mat.columns if c != "site_key"]
    if not probe_cols:
        raise ValueError("no probe columns")

    scored = mat.select(
        pl.col("site_key"),
        pl.sum_horizontal([pl.col(c) for c in probe_cols]).alias("total"),
        pl.max_horizontal([pl.col(c) for c in probe_cols]).alias("maxv"),
    ).with_columns((pl.col("maxv") / (pl.col("total") + 1e-9)).alias("specificity"))
    top = scored.sort("specificity", descending=True).head(top_n)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * top.height + 1)), constrained_layout=True)
    df = top.to_pandas()
    sns.barplot(data=df, y="site_key", x="specificity", color=preamble.COLOR_PRIDE, ax=ax, edgecolor="white")
    ax.invert_yaxis()
    ax.set_xlabel("Probe specificity (max / sum)")
    ax.set_ylabel("")
    ax.set_title("Chemoproteomic probe partitioning (PXD063995)")
    _finish_cat_ax(ax, horizontal=True)
    return _export(fig, "analysis_chemoproteomics")


def plot_triangulated_pxds(*, top_n: int = 15) -> dict[str, str]:
    tri = _load("triangulated_sites").head(top_n)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * tri.height + 1)), constrained_layout=True)
    df = tri.to_pandas()
    sns.barplot(
        data=df,
        y="site_key",
        x="n_pxds",
        color=preamble.COLOR_HIGHLIGHT,
        ax=ax,
        edgecolor="white",
    )
    ax.invert_yaxis()
    ax.set_xlabel("PRIDE deposits")
    ax.set_ylabel("")
    ax.set_title("Triangulated sites (canon + multi-PXD)")
    _finish_cat_ax(ax, horizontal=True)
    return _export(fig, "analysis_triangulated_pxds")


def plot_bap1ko_tissue() -> dict[str, str]:
    summary = _load("bap1ko_tissue_summary")
    df = summary.to_pandas()

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    sns.barplot(
        data=df,
        x="cond_tissue",
        y="n_sites",
        color=preamble.COLOR_CANON,
        edgecolor="white",
        ax=ax,
        width=0.6,
    )
    ax.set_ylabel("Unique O-GlcNAc sites")
    ax.set_xlabel("")
    ax.set_title("O-GlcNAc site burden by tissue (PXD035902)")
    _finish_cat_ax(ax)
    return _export(fig, "analysis_bap1ko_tissue")


def plot_gsea_contrast() -> dict[str, str]:
    panels: list[tuple[str, str, str]] = [
        ("canon_shared", "Canon-shared sites", preamble.COLOR_HIGHLIGHT),
        ("pride_novel", "PRIDE-novel sites", preamble.COLOR_PRIDE),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(_PANEL_W, 5), constrained_layout=True)
    any_panel = False
    for ax, (key, title, color) in zip(axes, panels, strict=True):
        path = ANALYSIS / f"pathway_{key}.parquet"
        if not path.is_file():
            ax.text(0.5, 0.5, "GSEA not run", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            continue
        any_panel = True
        pdf = pd.read_parquet(path).head(10)
        pdf["term"] = [str(t)[:45] + ("…" if len(str(t)) > 45 else "") for t in pdf["Term"]]
        pdf["score"] = -np.log10(pdf["Adjusted P-value"].clip(lower=1e-20))
        sns.barplot(data=pdf, y="term", x="score", color=color, edgecolor="white", ax=ax, zorder=3)
        ax.invert_yaxis()
        ax.set_ylabel("")
        ax.set_xlabel("−log₁₀(adj. P)")
        ax.set_title(title)
        _finish_cat_ax(ax, horizontal=True)
    if not any_panel:
        raise FileNotFoundError("run megatensor analyze with gseapy installed")
    fig.suptitle("Pathway enrichment contrasts (Enrichr GO BP)", fontsize=14, fontweight="bold")
    return _export(fig, "analysis_gsea_contrast")


def plot_silac_triangulation() -> dict[str, str]:
    fc = _load("silac_fc")
    tri_keys: set[str] = set()
    tri_path = ANALYSIS / "triangulated_sites.parquet"
    if tri_path.is_file():
        tri_keys = set(pl.read_parquet(tri_path)["site_key"].to_list())

    df = fc.to_pandas()
    df["triangulated"] = df["site_key"].isin(tri_keys)
    n_tri = int(df["triangulated"].sum())

    fig, ax = plt.subplots(figsize=(_PANEL_W, 6), constrained_layout=True)
    sns.scatterplot(
        data=df[~df["triangulated"]],
        x="log2_a",
        y="log2_fc",
        ax=ax,
        s=14,
        color=preamble.COLOR_MUTED,
        alpha=0.45,
        edgecolor=None,
        label="SILAC-only",
        rasterized=True,
    )
    sns.scatterplot(
        data=df[df["triangulated"]],
        x="log2_a",
        y="log2_fc",
        ax=ax,
        s=52,
        color=preamble.COLOR_HIGHLIGHT,
        edgecolor="white",
        linewidth=0.6,
        label=f"Triangulated ({n_tri})",
        zorder=4,
    )
    ax.axhline(0, color="#888888", lw=0.8, ls="--")
    ax.set_xlabel("log₂ A (mean abundance)")
    ax.set_ylabel("log₂(Heavy / Light)")
    ax.set_title(f"PXD039536 SILAC — {n_tri} triangulated sites among {fc.height} quantified")
    _finish_bar_ax(ax)
    ax.legend(loc="upper right", frameon=False)
    return _export(fig, "analysis_silac_triangulation")


def plot_silac_fc_histogram() -> dict[str, str]:
    fc = _load("silac_fc")
    log2 = fc["log2_fc"].to_numpy()
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    sns.histplot(log2, bins=30, ax=ax, color=preamble.COLOR_PRIDE, edgecolor="white", linewidth=0.6, alpha=0.85)
    med = float(np.median(log2))
    pct = 100 * (log2 > 0).mean()
    ax.axvline(med, color=preamble.COLOR_NEG, lw=1.4, ls="--", label=f"median {2**med:.2f}×")
    ax.axvline(0, color="#888888", lw=0.8, label="parity")
    ax.legend(title=f"{pct:.0f}% heavy-biased", frameon=False, loc="upper left")
    ax.set_xlabel("log₂(Heavy / Light)")
    ax.set_ylabel("Sites")
    ax.set_title("SILAC fold-change distribution (PXD039536)")
    _finish_bar_ax(ax)
    return _export(fig, "analysis_silac_fc_hist")


def plot_concordance_context() -> dict[str, str]:
    pairs = _load("concordance_annotated").filter(pl.col("n_shared_sites") >= 5)
    pdf = pairs.to_pandas()
    pdf["hue"] = pdf["comparison_class"].map(
        {
            "glycoid_family": "GlycoID family",
            "ogt_network": "OGT network",
            "ogt_cross_lab": "OGT cross-lab",
            "silac_cross_lab": "SILAC cross-lab",
            "cross_study": "Other",
        }
    ).fillna("Other")

    fig, ax = plt.subplots(figsize=(_PANEL_W, 5), constrained_layout=True)
    sns.barplot(
        data=pdf,
        y="comparison_label",
        x="pearson_r_log10",
        hue="hue",
        palette={
            "GlycoID family": preamble.COLOR_HIGHLIGHT,
            "OGT network": "#e07b39",
            "OGT cross-lab": preamble.COLOR_ACCENT,
            "SILAC cross-lab": preamble.COLOR_PRIDE,
            "Other": preamble.COLOR_MUTED,
        },
        dodge=False,
        ax=ax,
        edgecolor="white",
    )
    ax.invert_yaxis()
    ax.set_xlim(-0.05, 1.05)
    ax.set_xlabel("Pearson r")
    ax.set_ylabel("")
    ax.set_title("Within-family concordance exceeds cross-lab")
    ax.legend(title="Study context", frameon=False, loc="lower right")
    _finish_cat_ax(ax, horizontal=True)
    return _export(fig, "analysis_concordance_context")


def plot_ogt_concordance() -> dict[str, str]:
    pair = _load("concordance_ogt_pair")
    summary = json.loads((FIGURES / "analysis_summary.json").read_text())
    stats = summary.get("concordance_ogt", {})
    pdf = pair.to_pandas()

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    sns.scatterplot(data=pdf, x="log10_a", y="log10_b", ax=ax, s=32, color="#e07b39", alpha=0.7, edgecolor="white", linewidth=0.4)
    x, y = pdf["log10_a"].to_numpy(), pdf["log10_b"].to_numpy()
    if len(x) > 2:
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, coef[0] * xs + coef[1], color=preamble.COLOR_NEG, lw=1.2, ls="--")
    r = stats.get("pearson_r_log10", 0)
    n = stats.get("n_shared_sites", pair.height)
    ax.legend(
        handles=[Line2D([0], [0], color=preamble.COLOR_NEG, ls="--", label=f"r = {r}, n = {n}")],
        frameon=False,
        loc="lower right",
    )
    ax.set_xlabel("PXD035902 (BAP1KO glycomics)")
    ax.set_ylabel("PXD039536 (China SILAC OGT)")
    ax.set_title("OGT-network studies agree moderately at shared sites")
    _finish_bar_ax(ax)
    return _export(fig, "analysis_ogt_concordance")


def plot_bap1ko_brain_liver(*, label_top: int = 5) -> dict[str, str]:
    pairs = _load("bap1ko_brain_liver")
    pdf = pairs.to_pandas()
    x = pdf["log10_liver"].to_numpy()
    y = pdf["log10_brain"].to_numpy()

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    sns.scatterplot(data=pdf, x="log10_liver", y="log10_brain", ax=ax, s=16, color=preamble.COLOR_CANON, alpha=0.35, edgecolor=None, rasterized=True)
    lims = [min(x.min(), y.min()) - 0.2, max(x.max(), y.max()) + 0.2]
    ax.plot(lims, lims, color="#aaaaaa", lw=0.8, ls="--", zorder=1)
    order = np.argsort(pdf["brain_liver_ratio"].to_numpy())[::-1][:label_top]
    hi = pdf.iloc[order]
    sns.scatterplot(
        data=hi,
        x="log10_liver",
        y="log10_brain",
        ax=ax,
        s=48,
        color=preamble.COLOR_HIGHLIGHT,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )
    repel_labels(
        ax,
        hi["log10_liver"].to_numpy(),
        hi["log10_brain"].to_numpy(),
        hi["site_key"].tolist(),
        scatter_x=x,
        scatter_y=y,
        fontsize=7,
    )

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("log₁₀ liver intensity")
    ax.set_ylabel("log₁₀ brain intensity")
    ax.set_title(f"Brain vs liver ({pairs.height} shared sites)")
    _finish_bar_ax(ax)
    return _export(fig, "analysis_bap1ko_brain_liver")


def plot_evidence_ladder(*, top_n: int = 15) -> dict[str, str]:
    ladder = _load("site_evidence_ladder").head(top_n)
    pdf = ladder.to_pandas()

    fig, ax = plt.subplots(figsize=(8, max(4, 0.4 * len(pdf) + 1)), constrained_layout=True)
    sns.barplot(
        data=pdf,
        y="site_key",
        x="evidence_score",
        color=preamble.COLOR_HIGHLIGHT,
        edgecolor="white",
        ax=ax,
        zorder=3,
    )
    ax.invert_yaxis()
    ax.set_xlabel("Composite evidence score")
    ax.set_ylabel("")
    ax.set_title("Highest-evidence O-GlcNAc sites")
    _finish_cat_ax(ax, horizontal=True)
    return _export(fig, "analysis_evidence_ladder")


def plot_megatensor_impact() -> dict[str, str]:
    imp_path = FIGURES / "megatensor_importance.json"
    if not imp_path.is_file():
        imp_path = ANALYSIS / "megatensor_importance.json"
    if not imp_path.is_file():
        from megatensor.analysis.importance import write_importance_summary

        write_importance_summary()
    summary = json.loads(imp_path.read_text())
    payoff = summary["megatensor_payoff"]
    cost = summary["flat_file_cost"]
    pillars = summary["pillars"]
    n_pxd = pillars[0]["metrics"]["pxds_appended"]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)

    df_a = pd.DataFrame(
        {
            "metric": ["Pairwise PXD\nreconciliations", "Thin adapters\n(one per deposit)"],
            "count": [cost["pairwise_manual_joins"], n_pxd],
            "color": [preamble.COLOR_NEG, preamble.COLOR_HIGHLIGHT],
        }
    )
    sns.barplot(data=df_a, x="metric", y="count", hue="metric", palette=dict(zip(df_a["metric"], df_a["color"])), dodge=False, legend=False, ax=axes[0], edgecolor="white", width=0.55)
    axes[0].set_ylabel("Count")
    axes[0].set_xlabel("")
    axes[0].set_title("Flat files vs megatensor ingest")
    _finish_cat_ax(axes[0])

    df_b = pd.DataFrame(
        {
            "metric": ["Canon∩PRIDE\nsites", "Triangulated\nsites", "Multi-PXD\nsites"],
            "count": [payoff["shared_sites"], payoff["triangulated_sites"], pillars[2]["metrics"]["multi_pxd_sites"]],
        }
    )
    sns.barplot(data=df_b, x="metric", y="count", color=preamble.COLOR_CANON, edgecolor="white", ax=axes[1], width=0.55)
    axes[1].set_ylabel("Unique sites")
    axes[1].set_xlabel("")
    axes[1].set_title("Automatic cross-layer yield")
    _finish_cat_ax(axes[1])

    ctx = pillars[1]["metrics"]
    df_c = pd.DataFrame(
        {
            "metric": ["SILAC\nsites", "Tissue\nsites", "Chemo\nsites"],
            "count": [ctx["pride_sites_with_silac"], ctx["pride_sites_with_tissue"], ctx["chemoproteomic_sites"]],
        }
    )
    sns.barplot(data=df_c, x="metric", y="count", color=preamble.COLOR_PRIDE, edgecolor="white", ax=axes[2], width=0.55)
    axes[2].set_ylabel("Sites on condition axis")
    axes[2].set_xlabel("")
    axes[2].set_title("Biology stays sliceable")
    _finish_cat_ax(axes[2])

    fig.suptitle("Why the megatensor: structure unlocks cross-study O-GlcNAc science", fontsize=14, fontweight="bold")
    return _export(fig, "analysis_megatensor_impact")


def render_analysis_figures() -> dict[str, dict[str, str]]:
    FIGURES.mkdir(parents=True, exist_ok=True)
    builders = [
        ("megatensor_impact", plot_megatensor_impact),
        ("silac_ma", plot_silac_ma),
        ("silac_triangulation", plot_silac_triangulation),
        ("silac_fc_hist", plot_silac_fc_histogram),
        ("concordance_scatter", plot_concordance_scatter),
        ("concordance_context", plot_concordance_context),
        ("ogt_concordance", plot_ogt_concordance),
        ("replication_tiers", plot_replication_tiers),
        ("protein_hubs", plot_protein_hubs),
        ("evidence_ladder", plot_evidence_ladder),
        ("triangulated_pxds", plot_triangulated_pxds),
        ("chemoproteomics", plot_chemoproteomics),
        ("bap1ko_tissue", plot_bap1ko_tissue),
        ("bap1ko_brain_liver", plot_bap1ko_brain_liver),
        ("gsea_contrast", plot_gsea_contrast),
    ]
    paths: dict[str, dict[str, str]] = {}
    for key, fn in builders:
        try:
            paths[key] = fn()
        except Exception as exc:
            paths[key] = {"error": str(exc)}
    (FIGURES / "analysis_figures.json").write_text(json.dumps(paths, indent=2))
    return paths
