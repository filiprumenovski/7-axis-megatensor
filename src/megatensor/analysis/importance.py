"""Quantify why the megatensor representation matters (defense / preprint narrative)."""

from __future__ import annotations

import json
import math

import polars as pl

from megatensor.analysis.common import load_pride_obs, load_site_index
from megatensor.paths import FIGURES, ROOT
from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE


def megatensor_value_summary() -> dict:
    """Headline numbers tying architecture → capabilities → biology."""
    pride_summary = json.loads(PRIDE_STORE.summary_path.read_text()) if PRIDE_STORE.summary_path.is_file() else {}
    union_summary = json.loads(UNION_STORE.summary_path.read_text()) if UNION_STORE.summary_path.is_file() else {}
    canon_summary = json.loads(CANON_STORE.summary_path.read_text()) if CANON_STORE.summary_path.is_file() else {}

    datasets = pride_summary.get("datasets") or pride_summary.get("pxds") or []
    n_pxd = int(pride_summary["n_datasets"]) if "n_datasets" in pride_summary else len(datasets) or 12

    idx = load_site_index()
    obs = load_pride_obs()

    tri_path = ROOT / "megatensor" / "analysis" / "triangulated_sites.parquet"
    triangulated = pl.read_parquet(tri_path).height if tri_path.is_file() else 0

    multi_pxd = (
        pl.read_parquet(ROOT / "megatensor" / "analysis" / "multi_pxd_sites.parquet").height
        if (ROOT / "megatensor" / "analysis" / "multi_pxd_sites.parquet").is_file()
        else 0
    )

    shared = int(union_summary.get("shared_sites", idx.filter(pl.col("n_layers") >= 2).height))
    pride_sites = int(pride_summary.get("unique_sites", 0))
    canon_sites = int(canon_summary.get("unique_sites", 0))

    # Cross-axis coverage: same site_key usable across condition types
    site_key = pl.concat_str(
        [pl.col("protein_id_raw"), pl.col("residue_pos_raw").cast(pl.Utf8), pl.col("residue_aa")],
        separator=":",
    )
    obs = obs.with_columns(site_key.alias("site_key"))
    silac_sites = (
        obs.filter(pl.col("cond_treatment").is_in(["Light", "Heavy"]))
        .select("site_key")
        .unique()
        .height
    )
    tissue_sites = obs.filter(pl.col("cond_tissue").is_not_null()).select("site_key").unique().height
    chemo_sites = obs.filter(pl.col("dataset_id") == "PXD063995").select("site_key").unique().height

    silac_tri = 0
    silac_tri_path = ROOT / "megatensor" / "analysis" / "silac_fc.parquet"
    tri_keys_path = ROOT / "megatensor" / "analysis" / "triangulated_sites.parquet"
    if silac_tri_path.is_file() and tri_keys_path.is_file():
        silac_keys = set(pl.read_parquet(silac_tri_path)["site_key"].to_list())
        tri_keys = set(pl.read_parquet(tri_keys_path)["site_key"].to_list())
        silac_tri = len(silac_keys & tri_keys)

    pairwise_flat = n_pxd * (n_pxd - 1) // 2

    pillars = [
        {
            "pillar": "structural_interoperability",
            "claim": "One site identity joins canon and PRIDE without pairwise spreadsheet reconciliation.",
            "metrics": {
                "pxds_appended": n_pxd,
                "pairwise_reconciliations_avoided": pairwise_flat,
                "canon_pride_shared_sites": shared,
                "pct_pride_in_canon": round(100 * shared / max(pride_sites, 1), 1),
            },
        },
        {
            "pillar": "context_preservation",
            "claim": "Biology lives on condition axes; identity stays stable for cross-study slicing.",
            "metrics": {
                "pride_sites_with_silac": silac_sites,
                "pride_sites_with_tissue": tissue_sites,
                "chemoproteomic_sites": chemo_sites,
                "engines_represented": obs["source_engine"].drop_nulls().n_unique(),
                "countries_represented": obs["prov_country"].drop_nulls().n_unique(),
            },
        },
        {
            "pillar": "evidence_and_ml",
            "claim": "Replication tiers and exports turn deposits into ranked evidence and ML-ready tensors.",
            "metrics": {
                "multi_pxd_sites": multi_pxd,
                "triangulated_sites": triangulated,
                "silac_sites_also_triangulated": silac_tri,
                "observation_rows_after_set_rollup": obs.height,
                "ml_exports": ["exports/site_x_condition.parquet", "exports/site_x_features.parquet"],
            },
        },
    ]

    return {
        "purpose": "why_megatensor_matters",
        "headline": (
            "Separating site identity from experimental context lets 12 heterogeneous PRIDE deposits "
            "and two canon libraries coexist in one queryable object — enabling cross-study evidence "
            "ranking and condition-resolved biology that flat tables cannot express without N² reconciliation."
        ),
        "flat_file_cost": {
            "description": "Naïve flat-file meta-analysis requires harmonizing each PXD pair separately.",
            "pairwise_manual_joins": pairwise_flat,
            "canon_plus_pride_layers": 2,
        },
        "megatensor_payoff": {
            "adapters_only_at_ingest": True,
            "union_site_keys": idx.height,
            "canon_sites": canon_sites,
            "pride_sites": pride_sites,
            "shared_sites": shared,
            "triangulated_sites": triangulated,
        },
        "pillars": pillars,
        "example_queries_enabled": [
            "All observations for HCFC1 P51610:579:T across PXD035902, PXD039536, PXD058744",
            "SILAC Heavy/Light FC at sites that are also canon-supported and multi-PXD",
            "Brain vs liver intensity at the same site_key in BAP1KO glycomics",
            "Concordance at shared sites between China SILAC and US MaxQuant",
            "Chemoproteomic probe specificity matrix without re-parsing raw mzTab",
        ],
    }


def write_importance_summary() -> str:
    summary = megatensor_value_summary()
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / "megatensor_importance.json"
    path.write_text(json.dumps(summary, indent=2))
    analysis_root = ROOT / "megatensor" / "analysis"
    analysis_root.mkdir(parents=True, exist_ok=True)
    out = analysis_root / "megatensor_importance.json"
    out.write_text(json.dumps(summary, indent=2))
    return str(path)
