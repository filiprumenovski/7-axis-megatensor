"""Phase 4: identity enrichment (UniProt window/domain, optional disorder)."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.enrich.sasa import sasa_at_site
from megatensor.enrich.site_features import enrich_identities
from megatensor.paths import EXPORTS, FIGURES
from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE

log = structlog.get_logger()


def _union_identity_dim() -> pl.DataFrame:
    canon = pl.read_parquet(CANON_STORE.registry / "identity_dim.parquet")
    pride = pl.read_parquet(PRIDE_STORE.registry / "identity_dim.parquet")
    return pl.concat([canon, pride], how="diagonal_relaxed").unique(subset=["identity_id"])


def _run_gsea(genes: list[str], out_dir) -> dict | None:
    try:
        import gseapy as gp
    except ImportError:
        log.info("gseapy_skip", reason="not installed")
        return None

    genes = [g for g in genes if g]
    if len(genes) < 10:
        return None
    genes = list(dict.fromkeys(genes))[:500]
    try:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=["GO_Biological_Process_2023", "KEGG_2021_Human"],
            organism="human",
            outdir=str(out_dir / "gseapy"),
            no_plot=True,
        )
        if enr is not None and hasattr(enr, "results") and enr.results is not None:
            path = out_dir / "pathway_enrichment.parquet"
            enr.results.to_parquet(path)
            return {"genes": len(genes), "terms": int(enr.results.shape[0]), "path": str(path)}
    except Exception as exc:
        log.warning("gseapy_failed", error=str(exc))
    return None


def _figure_c_sasa() -> dict | None:
    traj_path = FIGURES / "figure_c_trajectory_candidate.csv"
    if not traj_path.is_file():
        return None
    import polars as pl

    traj = pl.read_csv(traj_path)
    if traj.is_empty():
        return None
    row = traj.row(0, named=True)
    acc = str(row["protein_id_raw"])
    pos = int(row["residue_pos_raw"])
    result = sasa_at_site(acc, pos)
    if result:
        (FIGURES / "figure_c_sasa.json").write_text(json.dumps(result, indent=2))
    return result


def run_enrich() -> dict:
    if not CANON_STORE.summary_path.exists() or not PRIDE_STORE.summary_path.exists():
        raise FileNotFoundError("run canon + pride-tensorize first")

    identity_dim = _union_identity_dim()
    enriched = enrich_identities(identity_dim)

    UNION_STORE.enrichment.mkdir(parents=True, exist_ok=True)
    out_path = UNION_STORE.enrichment / "site_features.parquet"
    enriched.write_parquet(out_path)

    FIGURES.mkdir(parents=True, exist_ok=True)
    completeness = {
        "sites": enriched.height,
        "seq_window_pct": round(100 * enriched["seq_match"].mean(), 2),
        "domain_pct": round(100 * enriched["region_type"].is_not_null().mean(), 2),
        "disorder_pct": round(100 * enriched["disorder_score"].is_not_null().mean(), 2),
        "gene_symbol_pct": round(100 * enriched["gene_symbol"].is_not_null().mean(), 2),
    }
    (FIGURES / "enrichment_completeness.json").write_text(json.dumps(completeness, indent=2))

    EXPORTS.mkdir(parents=True, exist_ok=True)
    gsea = _run_gsea(enriched["gene_symbol"].drop_nulls().to_list(), EXPORTS)

    sasa = _figure_c_sasa()

    summary = {
        "purpose": "identity_enrichment",
        "output": str(out_path),
        "completeness": completeness,
        "gsea": gsea,
        "figure_c_sasa": sasa,
    }
    (UNION_STORE.enrichment / "enrichment_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("enrichment_complete", **completeness)
    print(json.dumps(summary, indent=2))
    return summary
