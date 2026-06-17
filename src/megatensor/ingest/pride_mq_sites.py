"""MaxQuant O-GlcNAc / HexNAc(ST) site table adapter."""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from megatensor.ingest.pride_common import (
    base_provenance,
    collision_from_text,
    condition_from_name,
    finalize_obs,
    merge_conditions,
    parse_uniprot_acc,
    read_table,
)
from megatensor.ingest.pride_metrics import append_metric_rows, mq_metrics_for_row

INTENSITY_RE = re.compile(r"^Intensity", re.I)


def _pick_col(cols: list[str], *candidates: str) -> str | None:
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _condition_from_quant_col(col: str) -> dict[str, str | None]:
    return merge_conditions(condition_from_name(col))


def parse_mq_sites(path: Path, *, pxd: str, engine: str = "maxquant") -> pl.DataFrame:
    df = read_table(path)
    cols = df.columns
    pos_col = _pick_col(cols, "Position", "Positions within proteins", "Positions")
    prot_col = _pick_col(cols, "Protein", "Leading proteins", "Proteins")
    aa_col = _pick_col(cols, "Amino acid")
    loc_col = _pick_col(cols, "Localization prob", "Localization probability")
    raw_col = _pick_col(cols, "Best localization raw file", "Best score raw file")
    if not pos_col or not prot_col:
        return finalize_obs(pl.DataFrame())

    file_cond = condition_from_name(path.name)
    prov = base_provenance(pxd, engine, path)
    intensity_cols = [c for c in cols if INTENSITY_RE.match(c)]

    rows: list[dict] = []
    for row in df.iter_rows(named=True):
        acc = parse_uniprot_acc(str(row.get(prot_col) or ""))
        pos = row.get(pos_col)
        try:
            pos_i = int(float(pos)) if pos is not None and str(pos) != "" else None
        except (TypeError, ValueError):
            pos_i = None
        aa = str(row.get(aa_col) or "").strip().upper()[:1] if aa_col else None
        if not acc or not pos_i or aa not in ("S", "T"):
            continue
        loc = row.get(loc_col) if loc_col else None
        try:
            loc_f = float(loc) if loc not in (None, "") else None
        except (TypeError, ValueError):
            loc_f = None

        raw_file = str(row.get(raw_col) or "") if raw_col else ""
        row_cond = merge_conditions(file_cond, condition_from_name(raw_file))
        collision = collision_from_text(raw_file)

        base = {
            **prov,
            **row_cond,
            "protein_id_raw": acc,
            "protein_id_type": "uniprot_acc",
            "isoform_raw": None,
            "residue_pos_raw": pos_i,
            "residue_aa": aa,
            "loc_score": loc_f,
            "loc_method": "mq_locprob",
            "loc_is_ambiguous": bool(loc_f is not None and loc_f < 0.75),
            "acq_ms_mode": "DDA",
            "acq_msn_level": "MS2",
            "acq_gradient_min": None,
            "acq_collision": collision,
            "acq_enrichment": None,
            "cond_replicate": raw_file or None,
            "qc_flags": None,
        }

        if intensity_cols:
            for ic in intensity_cols:
                col_cond = merge_conditions(row_cond, _condition_from_quant_col(ic))
                metrics = mq_metrics_for_row(row, ic)
                if not metrics:
                    continue
                append_metric_rows(rows, {**base, **col_cond, "cond_replicate": ic}, metrics)
        else:
            append_metric_rows(rows, base, [("spectral_count", 1.0, "curated", "MQ_site")])

    return finalize_obs(pl.DataFrame(rows))
