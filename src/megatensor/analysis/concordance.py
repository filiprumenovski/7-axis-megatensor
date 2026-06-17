"""Cross-PXD intensity concordance at shared sites."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import load_pride_obs


def intensity_by_site_pxd() -> pl.DataFrame:
    return (
        load_pride_obs()
        .filter(pl.col("metric_name") == "intensity")
        .with_columns(
            pl.concat_str(
                [pl.col("protein_id_raw"), pl.col("residue_pos_raw").cast(pl.Utf8), pl.col("residue_aa")],
                separator=":",
            ).alias("site_key")
        )
        .group_by(["site_key", "dataset_id"])
        .agg(pl.col("metric_value").mean().alias("intensity_mean"), pl.len().alias("n_obs"))
    )


def pairwise_concordance(
    pxd_a: str,
    pxd_b: str,
    *,
    min_obs: int = 1,
) -> tuple[pl.DataFrame, dict]:
    wide = intensity_by_site_pxd()
    a = wide.filter(pl.col("dataset_id") == pxd_a).rename({"intensity_mean": "a", "n_obs": "n_a"})
    b = wide.filter(pl.col("dataset_id") == pxd_b).rename({"intensity_mean": "b", "n_obs": "n_b"})
    joined = a.join(b, on="site_key", how="inner").filter(pl.col("n_a") >= min_obs, pl.col("n_b") >= min_obs)
    if joined.is_empty():
        return joined, {}

    import numpy as np

    la = np.log10(joined["a"].to_numpy() + 1.0)
    lb = np.log10(joined["b"].to_numpy() + 1.0)
    r = float(np.corrcoef(la, lb)[0, 1]) if len(la) > 1 else 0.0
    rho = r
    try:
        from scipy.stats import spearmanr

        rho = float(spearmanr(la, lb).statistic)
    except Exception:
        pass

    stats = {
        "pxd_a": pxd_a,
        "pxd_b": pxd_b,
        "n_shared_sites": joined.height,
        "pearson_r_log10": round(r, 3),
        "spearman_rho_log10": round(rho, 3),
    }
    return joined.with_columns(pl.Series("log10_a", la), pl.Series("log10_b", lb)), stats


def all_pairwise_concordance(pxds: list[str] | None = None) -> pl.DataFrame:
    wide = intensity_by_site_pxd()
    pxds = pxds or wide["dataset_id"].unique().sort().to_list()
    rows: list[dict] = []
    for i, a in enumerate(pxds):
        for b in pxds[i + 1 :]:
            _, stats = pairwise_concordance(a, b)
            if stats:
                rows.append(stats)
    return pl.DataFrame(rows).sort("n_shared_sites", descending=True)


# Narrative labels for concordance pairs (engine / study family).
_PAIR_CONTEXT: dict[frozenset[str], tuple[str, str]] = {
    frozenset({"PXD033026", "PXD033043"}): ("glycoid_family", "GlycoID cytosol × insulin serum"),
    frozenset({"PXD033043", "PXD033062"}): ("glycoid_family", "GlycoID insulin × proximity serum"),
    frozenset({"PXD033026", "PXD033062"}): ("glycoid_family", "GlycoID cytosol × proximity serum"),
    frozenset({"PXD035902", "PXD039536"}): ("ogt_network", "OGT interactome × SILAC OGT"),
    frozenset({"PXD035902", "PXD058744"}): ("ogt_cross_lab", "BAP1KO glycomics × US MaxQuant"),
    frozenset({"PXD039536", "PXD058744"}): ("silac_cross_lab", "China SILAC × US MaxQuant"),
}


def annotate_concordance_pairs(pairs: pl.DataFrame) -> pl.DataFrame:
    """Add comparison class and short label to pairwise concordance rows."""
    if pairs.is_empty():
        return pairs

    classes: list[str] = []
    labels: list[str] = []
    for r in pairs.iter_rows(named=True):
        key = frozenset({r["pxd_a"], r["pxd_b"]})
        cls, lbl = _PAIR_CONTEXT.get(key, ("cross_study", f"{r['pxd_a']} × {r['pxd_b']}"))
        classes.append(cls)
        labels.append(lbl)
    return pairs.with_columns(pl.Series("comparison_class", classes), pl.Series("comparison_label", labels))


def concordance_context_summary(pairs: pl.DataFrame) -> dict:
    """Aggregate r by comparison class (within-pipeline vs cross-lab)."""
    ann = annotate_concordance_pairs(pairs)
    if ann.is_empty():
        return {}
    by_class = (
        ann.group_by("comparison_class")
        .agg(
            pl.len().alias("n_pairs"),
            pl.col("pearson_r_log10").mean().alias("mean_r"),
            pl.col("n_shared_sites").sum().alias("total_shared_sites"),
        )
        .sort("mean_r", descending=True)
    )
    return {
        "pairs_annotated": ann.height,
        "by_class": by_class.to_dicts(),
        "glycoid_replicate_r": ann.filter(pl.col("comparison_class") == "glycoid_family")["pearson_r_log10"].mean(),
        "cross_lab_mean_r": ann.filter(pl.col("comparison_class").is_in(["silac_cross_lab", "ogt_cross_lab"]))[
            "pearson_r_log10"
        ].mean(),
    }
