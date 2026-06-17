"""Proteome Discoverer flattened Excel export (peptide + site subtables)."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from megatensor.ingest.pride_common import (
    base_provenance,
    condition_from_name,
    finalize_obs,
    merge_conditions,
    parse_uniprot_acc,
)
from megatensor.ingest.pride_metrics import append_metric_rows


def _read_flat_sheet(path: Path) -> pl.DataFrame:
    return pl.read_excel(path, has_header=False)


def parse_pd_xlsx(path: Path, *, pxd: str) -> pl.DataFrame:
    df = _read_flat_sheet(path)
    file_cond = condition_from_name(path.name)
    prov = base_provenance(pxd, "pd", path)
    rows: list[dict] = []

    for vals in df.iter_rows():
        cells = [str(v).strip() if v is not None else "" for v in vals]
        if len(cells) < 15:
            continue
        if cells[4] != "HexNAc":
            continue
        aa = cells[5].upper()[:1]
        if aa not in ("S", "T"):
            continue
        try:
            loc = float(cells[7]) if cells[7] else None
            pos = int(float(cells[14])) if cells[14] else None
        except (TypeError, ValueError):
            continue
        acc = parse_uniprot_acc(cells[12])
        if not acc or not pos:
            continue
        if loc is not None and loc > 1:
            loc = loc / 100.0

        base = {
            **prov,
            **merge_conditions(file_cond),
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
            "acq_collision": None,
            "acq_enrichment": None,
            "cond_replicate": None,
            "qc_flags": None,
        }
        metrics = [("score", loc, "raw", "ptmRS")] if loc is not None else [("spectral_count", 1.0, "curated", "PD_site")]
        append_metric_rows(rows, base, metrics)

    return finalize_obs(pl.DataFrame(rows))
