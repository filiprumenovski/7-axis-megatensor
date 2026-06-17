"""Orchestrate cross-layer analyses → tables + summary JSON."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.analysis.bap1ko import bap1ko_tissue_summary, bap1ko_top_sites
from megatensor.analysis.conditions import bap1ko_tissue_sites, chemoproteomics_site_matrix
from megatensor.analysis.concordance import (
    all_pairwise_concordance,
    annotate_concordance_pairs,
    concordance_context_summary,
    pairwise_concordance,
)
from megatensor.analysis.evidence import site_evidence_ladder, triangulated_intensity_matrix
from megatensor.analysis.hubs import hubs_with_genes
from megatensor.analysis.importance import megatensor_value_summary, write_importance_summary
from megatensor.analysis.novel import novel_summary, site_gene_lists
from megatensor.analysis.replication import replication_tables
from megatensor.analysis.silac import silac_fold_changes, silac_summary
from megatensor.analysis.tissue_contrast import bap1ko_brain_liver_pairs, bap1ko_brain_liver_summary
from megatensor.paths import FIGURES, ROOT
from megatensor.pipeline.supplementary import write_supplementary
from megatensor.store import UNION_STORE

log = structlog.get_logger()

ANALYSIS_ROOT = ROOT / "megatensor" / "analysis"


def _write(df: pl.DataFrame, name: str) -> str:
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)
    path = ANALYSIS_ROOT / f"{name}.parquet"
    if df.is_empty():
        return ""
    df.write_parquet(path)
    return str(path)


def _run_gsea(gene_lists: dict[str, list[str]]) -> dict | None:
    try:
        import gseapy as gp
    except ImportError:
        log.info("gseapy_skip", reason="pip install gseapy")
        return None

    out: dict = {}
    for label, genes in gene_lists.items():
        genes = list(dict.fromkeys(g for g in genes if g))[:800]
        if len(genes) < 15:
            log.info("gsea_skip_small", label=label, n=len(genes))
            continue
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human"],
                organism="human",
                outdir=str(ANALYSIS_ROOT / f"gseapy_{label}"),
                no_plot=True,
            )
            if enr is not None and enr.results is not None and not enr.results.empty:
                full_path = ANALYSIS_ROOT / f"pathway_{label}_full.parquet"
                enr.results.to_parquet(full_path)
                top = enr.results.head(20)
                top_path = ANALYSIS_ROOT / f"pathway_{label}.parquet"
                top.to_parquet(top_path)
                out[label] = {
                    "genes": len(genes),
                    "terms": int(enr.results.shape[0]),
                    "path": str(top_path),
                    "top_term": str(enr.results.iloc[0].get("Term", "")),
                    "top_padj": float(enr.results.iloc[0].get("Adjusted P-value", 1)),
                }
        except Exception as exc:
            log.warning("gsea_failed", label=label, error=str(exc))
    return out or None


def run_gsea_step() -> dict | None:
    """Pathway enrichment after `enrich` has written site_features."""
    genes = site_gene_lists()
    gsea = _run_gsea(genes) if genes else None
    summary_path = ANALYSIS_ROOT / "analysis_summary.json"
    fig_path = FIGURES / "analysis_summary.json"
    if summary_path.is_file():
        summary = json.loads(summary_path.read_text())
        summary["gsea"] = gsea
        summary_path.write_text(json.dumps(summary, indent=2, default=str))
        fig_path.write_text(json.dumps(summary, indent=2, default=str))
    return gsea


def run_analyze(*, skip_gsea: bool = False) -> dict:
    if not (UNION_STORE.staging / "site_index.parquet").is_file():
        raise FileNotFoundError("run union first")

    rep = replication_tables()
    hubs = hubs_with_genes(rep["protein_hubs"])
    silac = silac_fold_changes()
    concordance_pair, conc_stats = pairwise_concordance("PXD039536", "PXD058744")
    concordance_all = all_pairwise_concordance()
    concordance_annotated = annotate_concordance_pairs(concordance_all)
    concordance_context = concordance_context_summary(concordance_all)
    ogt_pair, ogt_stats = pairwise_concordance("PXD035902", "PXD039536")
    evidence = site_evidence_ladder(top_n=50)
    tri_heatmap = triangulated_intensity_matrix(top_n=25)
    brain_liver = bap1ko_brain_liver_pairs()
    brain_liver_stats = bap1ko_brain_liver_summary(brain_liver)
    importance = megatensor_value_summary()
    importance_path = write_importance_summary()
    novel = novel_summary()
    chemo = chemoproteomics_site_matrix()
    bap1 = bap1ko_tissue_sites()
    bap1_summary = bap1ko_tissue_summary()
    bap1_top = bap1ko_top_sites()
    genes = site_gene_lists() if not skip_gsea else {}
    gsea = _run_gsea(genes) if genes and not skip_gsea else None
    supp = write_supplementary()

    paths = {
        "multi_pxd_sites": _write(rep["multi_pxd_sites"], "multi_pxd_sites"),
        "triangulated_sites": _write(rep["triangulated_sites"], "triangulated_sites"),
        "protein_hubs": _write(hubs, "protein_hubs"),
        "overlap_tiers": _write(rep["overlap_tiers"], "overlap_tiers"),
        "silac_fc": _write(silac, "silac_fc"),
        "concordance_silac_pair": _write(concordance_pair, "concordance_silac_pair"),
        "concordance_all_pairs": _write(concordance_all, "concordance_all_pairs"),
        "concordance_annotated": _write(concordance_annotated, "concordance_annotated"),
        "concordance_ogt_pair": _write(ogt_pair, "concordance_ogt_pair"),
        "site_evidence_ladder": _write(evidence, "site_evidence_ladder"),
        "triangulated_intensity_matrix": _write(tri_heatmap, "triangulated_intensity_matrix"),
        "bap1ko_brain_liver": _write(brain_liver, "bap1ko_brain_liver"),
        "novel_summary": _write(novel, "novel_summary"),
        "chemoproteomics_matrix": _write(chemo, "chemoproteomics_matrix"),
        "bap1ko_tissue": _write(bap1, "bap1ko_tissue"),
        "bap1ko_tissue_summary": _write(bap1_summary, "bap1ko_tissue_summary"),
        "bap1ko_top_sites": _write(bap1_top, "bap1ko_top_sites"),
    }

    top_hub = hubs.row(0, named=True) if hubs.height else None
    summary = {
        "purpose": "cross_layer_analysis",
        "replication": {
            "multi_pxd_sites": rep["multi_pxd_sites"].height,
            "triangulated_sites": rep["triangulated_sites"].height,
            "top_protein_hub": top_hub,
        },
        "silac": silac_summary(silac),
        "concordance": conc_stats,
        "concordance_context": concordance_context,
        "concordance_ogt": ogt_stats,
        "brain_liver": brain_liver_stats,
        "importance": importance,
        "concordance_pairs": concordance_annotated.head(10).to_dicts() if concordance_annotated.height else [],
        "novel": novel.to_dicts(),
        "gsea": gsea,
        "supplementary": supp,
        "outputs": {k: v for k, v in paths.items() if v},
        "importance_path": importance_path,
    }

    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=True)
    (ANALYSIS_ROOT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    FIGURES.mkdir(parents=True, exist_ok=True)
    (FIGURES / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    log.info("analysis_complete", **{k: v for k, v in summary.items() if k not in ("outputs", "concordance_pairs")})
    print(json.dumps(summary, indent=2, default=str))
    return summary
