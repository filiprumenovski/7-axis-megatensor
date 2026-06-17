"""Build per-site enrichment table from UniProt (+ optional disorder)."""

from __future__ import annotations

import polars as pl

from megatensor.enrich.uniprot import (
    _normalize_acc,
    features_at_position,
    fetch_batch,
    gene_symbol,
    sequence_of,
)

WINDOW = 7


def _seq_window(seq: str, pos: int, aa: str, k: int = WINDOW) -> str | None:
    if not seq or pos < 1 or pos > len(seq):
        return None
    if seq[pos - 1].upper() != aa.upper():
        return None
    lo = max(0, pos - 1 - k)
    hi = min(len(seq), pos + k)
    return seq[lo:hi]


def _region_label(features: list[dict]) -> tuple[str | None, str | None]:
    if not features:
        return None, None
    priority = ("Domain", "Repeat", "Region", "Motif", "Chain")
    for p in priority:
        for f in features:
            if f.get("type") == p:
                return p, f.get("description")
    f0 = features[0]
    return f0.get("type"), f0.get("description")


def _disorder_scores(sequences: dict[str, str]) -> dict[str, list[float]]:
    try:
        import metapredict
    except ImportError:
        return {}

    out: dict[str, list[float]] = {}
    for acc, seq in sequences.items():
        if not seq:
            continue
        try:
            scores = metapredict.predict_disorder(seq)
            if scores is not None and len(scores) == len(seq):
                out[acc] = [float(x) for x in scores]
        except Exception:
            continue
    return out


def enrich_identities(identity_dim: pl.DataFrame) -> pl.DataFrame:
    """Attach seq_window, domain, gene_symbol, optional disorder_score per site."""
    sites = identity_dim.filter(~pl.col("protein_level_only") & pl.col("protein_acc").is_not_null())
    accs = [_normalize_acc(a) for a in sites["protein_acc"].unique().sort().to_list()]
    accs = [a for a in accs if a]
    entries = fetch_batch(accs)

    sequences = {acc: sequence_of(e) for acc, e in entries.items()}
    disorder = _disorder_scores(sequences)

    rows: list[dict] = []
    for row in sites.iter_rows(named=True):
        acc = _normalize_acc(row["protein_acc"])
        pos_raw = row["residue_pos"]
        aa = row["residue_aa"]
        if acc is None or pos_raw is None or aa not in ("S", "T"):
            continue
        pos = int(pos_raw)
        entry = entries.get(acc, {})
        seq = sequences.get(acc, "")
        feats = features_at_position(entry, pos) if entry else []
        region_type, domain_name = _region_label(feats)
        window = _seq_window(seq, pos, aa)
        d_score = None
        if acc in disorder and 0 < pos <= len(disorder[acc]):
            d_score = disorder[acc][pos - 1]

        rows.append(
            {
                "identity_id": row["identity_id"],
                "protein_acc": acc,
                "residue_pos": pos,
                "residue_aa": aa,
                "gene_symbol": gene_symbol(entry) if entry else None,
                "seq_window": window,
                "region_type": region_type,
                "domain_name": domain_name,
                "disorder_score": d_score,
                "seq_match": window is not None,
            }
        )

    return pl.DataFrame(rows)
