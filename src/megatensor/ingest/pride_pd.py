"""Proteome Discoverer PSM export adapter (Byonic / Sequest)."""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from megatensor.ingest.pride_common import (
    GLYCO_MOD_RE,
    SITE_MOD_RE,
    base_provenance,
    collision_from_text,
    condition_from_name,
    finalize_obs,
    merge_conditions,
    parse_uniprot_acc,
    read_table,
)
from megatensor.ingest.pride_metrics import append_metric_rows, pd_metrics_for_row

PTMRS_SITE_RE = re.compile(r"([ST])\((\d+)\):\s*([0-9.]+)")


def _parse_mod_sites(mod: str, pos_in_protein: int | None) -> list[tuple[str, int]]:
    sites: list[tuple[str, int]] = []
    if not mod or not GLYCO_MOD_RE.search(mod):
        return sites
    glyco_hits = [m for m in SITE_MOD_RE.finditer(mod) if GLYCO_MOD_RE.search(m.group(3))]
    if pos_in_protein is not None and glyco_hits:
        sites.append((glyco_hits[0].group(1), pos_in_protein))
    return sites


def _best_loc_score(row: dict, cols: dict[str, str]) -> float | None:
    scores: list[float] = []
    for key, col in cols.items():
        if "ptmrs" not in key:
            continue
        if not any(tok in key for tok in ("site", "glcnac", "hexnac", "best")):
            continue
        text = str(row.get(col) or "")
        for m in PTMRS_SITE_RE.finditer(text):
            try:
                scores.append(float(m.group(3)))
            except ValueError:
                pass
    return max(scores) if scores else None


def _infer_search_node(cols: dict[str, str], row: dict) -> str | None:
    node = row.get(cols.get("identifying node", "")) or row.get(cols.get("search id", ""))
    return str(node)[:64] if node else None


def parse_pd_psm(path: Path, *, pxd: str) -> pl.DataFrame:
    df = read_table(path)
    cols = {c.lower(): c for c in df.columns}
    mod_col = cols.get("modifications")
    acc_col = cols.get("master protein accessions") or cols.get("protein accessions")
    pos_col = cols.get("position in protein") or cols.get("positions in master proteins")
    activation_col = cols.get("activation type")
    if not mod_col or not acc_col:
        return finalize_obs(pl.DataFrame())

    file_cond = condition_from_name(path.name)
    prov = base_provenance(pxd, "pd", path)
    rows: list[dict] = []

    for row in df.iter_rows(named=True):
        mod = str(row.get(mod_col) or "")
        if not GLYCO_MOD_RE.search(mod):
            continue
        acc = parse_uniprot_acc(str(row.get(acc_col) or ""))
        if not acc:
            continue
        pos_raw = row.get(pos_col) if pos_col else None
        try:
            pos_in_protein = int(float(pos_raw)) if pos_raw not in (None, "") else None
        except (TypeError, ValueError):
            pos_in_protein = None

        spectrum = str(row.get(cols.get("spectrum file", "")) or "")
        if not prov.get("prov_search_version"):
            node = _infer_search_node(cols, row)
            if node:
                prov = {**prov, "prov_search_version": node}
        row_cond = merge_conditions(
            file_cond,
            condition_from_name(spectrum),
        )
        collision = collision_from_text(str(row.get(activation_col) or "")) or collision_from_text(spectrum)

        for aa, pos in _parse_mod_sites(mod, pos_in_protein):
            if pos_in_protein is None:
                continue
            loc = _best_loc_score(row, cols)
            metrics = pd_metrics_for_row(row, cols)
            base = {
                **prov,
                **row_cond,
                "protein_id_raw": acc,
                "protein_id_type": "uniprot_acc",
                "isoform_raw": None,
                "residue_pos_raw": pos,
                "residue_aa": aa,
                "loc_score": loc,
                "loc_method": "ptmrs",
                "loc_is_ambiguous": bool(loc is not None and loc < 0.75),
                "acq_ms_mode": "DDA",
                "acq_msn_level": "MS2",
                "acq_gradient_min": None,
                "acq_collision": collision,
                "acq_enrichment": None,
                "cond_replicate": None,
                "qc_flags": None,
            }
            append_metric_rows(rows, base, metrics)

    return finalize_obs(pl.DataFrame(rows))
