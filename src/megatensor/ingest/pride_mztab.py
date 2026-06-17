"""mzTab (PD export) PSM adapter."""

from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from megatensor.ingest.pride_common import (
    GLYCO_MOD_RE,
    base_provenance,
    collision_from_text,
    condition_from_name,
    finalize_obs,
    merge_conditions,
    parse_uniprot_acc,
)
from megatensor.ingest.pride_metrics import append_metric_rows, mztab_metrics_for_row

MOD_RE = re.compile(r"(\d+)-UNIMOD:(\d+)")
GLYCO_UNIMOD = {"43", "934", "121"}


def _parse_mztab_psm(path: Path) -> list[dict]:
    header: list[str] | None = None
    rows: list[dict] = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line or line.startswith("MTD") or line.startswith("COM"):
            continue
        parts = line.split("\t")
        if parts[0] == "PSH":
            header = parts[1:]
            continue
        if parts[0] != "PSM" or not header:
            continue
        rows.append(dict(zip(header, parts[1:], strict=False)))
    return rows


def parse_mztab(path: Path, *, pxd: str) -> pl.DataFrame:
    prov = base_provenance(pxd, "pd", path)
    file_cond = condition_from_name(path.name)
    for label in ("PC", "DDE", "DADPS", "DIAZO"):
        if label.lower() in path.stem.lower():
            file_cond = merge_conditions(file_cond, {"cond_treatment": label})
            break

    out_rows: list[dict] = []
    for row in _parse_mztab_psm(path):
        mods = row.get("modifications") or ""
        if not GLYCO_MOD_RE.search(mods):
            continue
        acc = parse_uniprot_acc(row.get("accession"))
        if not acc:
            continue
        try:
            start = int(float(row.get("start") or 0))
        except (TypeError, ValueError):
            continue
        seq = row.get("sequence") or ""
        spectra = str(row.get("spectra_ref") or "")
        row_cond = merge_conditions(file_cond, condition_from_name(spectra))
        collision = collision_from_text(path.name) or "HCD"

        for m in MOD_RE.finditer(mods):
            unimod = m.group(2)
            if unimod not in GLYCO_UNIMOD and not GLYCO_MOD_RE.search(mods):
                continue
            pep_pos = int(m.group(1))
            if pep_pos < 1 or pep_pos > len(seq):
                continue
            aa = seq[pep_pos - 1].upper()
            if aa not in ("S", "T"):
                continue
            protein_pos = start + pep_pos - 1
            try:
                loc = float(row.get("search_engine_score[2]") or row.get("search_engine_score[1]") or 0)
            except (TypeError, ValueError):
                loc = None
            base = {
                **prov,
                **row_cond,
                "protein_id_raw": acc,
                "protein_id_type": "uniprot_acc",
                "isoform_raw": None,
                "residue_pos_raw": protein_pos,
                "residue_aa": aa,
                "loc_score": loc,
                "loc_method": "ptmrs",
                "loc_is_ambiguous": bool(loc is not None and loc < 0.75),
                "acq_ms_mode": "DDA",
                "acq_msn_level": "MS2",
                "acq_gradient_min": None,
                "acq_collision": collision,
                "acq_enrichment": "chemoenzymatic",
                "cond_replicate": spectra or None,
                "qc_flags": None,
            }
            append_metric_rows(out_rows, base, mztab_metrics_for_row(row))

    return finalize_obs(pl.DataFrame(out_rows))
