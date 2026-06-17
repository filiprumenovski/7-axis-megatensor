"""Build per-site enrichment table for identities in the megatensor."""

from __future__ import annotations

import polars as pl
import structlog

from megatensor.enrich.disorder import disorder_at_site
from megatensor.enrich.uniprot import feature_at_site, fetch_accessions, gene_symbol, seq_window, sequence_for
from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE

log = structlog.get_logger()


def load_site_identities() -> pl.DataFrame:
    frames = []
    for store in (CANON_STORE, PRIDE_STORE):
        path = store.registry / "identity_dim.parquet"
        if path.is_file():
            frames.append(pl.read_parquet(path).filter(~pl.col("protein_level_only")))
    if not frames:
        raise FileNotFoundError("No identity_dim — run canon + pride-tensorize first")
    return pl.concat(frames, how="diagonal_relaxed").unique("identity_id")


def enrich_sites(*, flank: int = 7) -> pl.DataFrame:
    sites = load_site_identities()
    accs = sites["protein_acc"].drop_nulls().unique().sort().to_list()
    log.info("enrich_uniprot_fetch", proteins=len(accs))
    entries = fetch_accessions(accs)

    rows: list[dict] = []
    for site in sites.iter_rows(named=True):
        acc = site["protein_acc"]
        pos = site["residue_pos"]
        entry = entries.get(acc) if acc else None
        seq = sequence_for(entry) if entry else None
        region_type, domain_name = feature_at_site(entry, int(pos)) if entry and pos else (None, None)
        disorder = disorder_at_site(seq, int(pos)) if seq and pos else None
        window = seq_window(seq, int(pos), flank=flank) if seq and pos else None
        rows.append(
            {
                "identity_id": site["identity_id"],
                "protein_acc": acc,
                "residue_pos": pos,
                "residue_aa": site["residue_aa"],
                "gene_symbol": gene_symbol(entry) if entry else None,
                "seq_window": window,
                "region_type": region_type,
                "domain_name": domain_name,
                "disorder_score": disorder,
                "disorder_region": bool(disorder is not None and disorder >= 0.5),
                "uniprot_resolved": entry is not None,
            }
        )

    out = pl.DataFrame(rows)
    UNION_STORE.enrichment.mkdir(parents=True, exist_ok=True)
    out.write_parquet(UNION_STORE.enrichment / "site_features.parquet")
    summary = {
        "sites": out.height,
        "uniprot_resolved": int(out["uniprot_resolved"].sum()),
        "with_window": int(out["seq_window"].is_not_null().sum()),
        "with_disorder": int(out["disorder_score"].is_not_null().sum()),
        "with_domain": int(out["region_type"].is_not_null().sum()),
    }
    log.info("enrich_sites_complete", **summary)
    return out
