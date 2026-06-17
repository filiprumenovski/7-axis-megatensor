"""DIA-NN report.tsv adapter (site-localized O-GlcNAc precursors)."""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import polars as pl

from megatensor.ingest.pride_common import (
    GLYCO_MOD_RE,
    base_provenance,
    finalize_obs,
    parse_uniprot_acc,
)

# DIA-NN encodes mods as S(UniMod:43) or T(UniMod:934)
DIANN_MOD_RE = re.compile(r"([ST])\(UniMod:(43|934|121)\)", re.I)


def _report_path_in_zip(path: Path) -> str:
    with zipfile.ZipFile(path) as zf:
        names = [n for n in zf.namelist() if n.endswith("report.tsv") and "first-pass" not in n]
        if not names:
            names = [n for n in zf.namelist() if n.endswith("report.tsv")]
        if not names:
            raise FileNotFoundError(f"No report.tsv in {path}")
        return sorted(names, key=len)[0]


def parse_diann_zip(path: Path, *, pxd: str) -> pl.DataFrame:
    prov = base_provenance(pxd, "diann", path)
    member = _report_path_in_zip(path)
    rows: list[dict] = []

    with zipfile.ZipFile(path) as zf, zf.open(member) as fh:
        header = fh.readline().decode(errors="ignore").rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        need = ["Modified.Sequence", "Protein.Ids", "Run", "Precursor.Quantity", "PTM.Site.Confidence"]
        if not all(c in idx for c in need[:3]):
            return finalize_obs(pl.DataFrame())
        for raw in fh:
            parts = raw.decode(errors="ignore").rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            modseq = parts[idx["Modified.Sequence"]]
            if not DIANN_MOD_RE.search(modseq):
                continue
            acc = parse_uniprot_acc(parts[idx["Protein.Ids"]].split(";")[0])
            if not acc:
                continue
            stripped = parts[idx.get("Stripped.Sequence", idx["Modified.Sequence"])]
            for m in DIANN_MOD_RE.finditer(modseq):
                aa = m.group(1).upper()
                # site index in stripped peptide (1-based)
                pos = _mod_position(stripped, modseq, m.start())
                if pos is None:
                    continue
                try:
                    qty = float(parts[idx["Precursor.Quantity"]])
                except (TypeError, ValueError):
                    qty = 1.0
                try:
                    loc = float(parts[idx.get("PTM.Site.Confidence", idx["Modified.Sequence"])])
                except (TypeError, ValueError, IndexError):
                    loc = None
                rows.append(
                    {
                        **prov,
                        "protein_id_raw": acc,
                        "protein_id_type": "uniprot_acc",
                        "isoform_raw": None,
                        "residue_pos_raw": pos,
                        "residue_aa": aa,
                        "loc_score": loc,
                        "loc_method": "diann_ptm",
                        "loc_is_ambiguous": bool(loc is not None and loc < 0.75),
                        "cond_cell_line": None,
                        "cond_tissue": None,
                        "cond_treatment": None,
                        "cond_timepoint": None,
                        "cond_fraction": None,
                        "cond_replicate": parts[idx["Run"]],
                        "cond_replicate_type": None,
                        "acq_ms_mode": "DIA",
                        "acq_msn_level": "MS2",
                        "acq_gradient_min": None,
                        "acq_collision": None,
                        "acq_enrichment": None,
                        "metric_name": "intensity",
                        "metric_value": qty,
                        "metric_norm_state": "raw",
                        "metric_unit": None,
                        "qc_flags": "peptide_position_only",
                    }
                )

    return finalize_obs(pl.DataFrame(rows))


def _mod_position(stripped: str, modified: str, mod_start: int) -> int | None:
    """Map mod token offset in Modified.Sequence to 1-based index in stripped peptide."""
    # count AA characters before mod_start in modified string (skip parenthetical mods)
    aa_idx = 0
    i = 0
    while i < mod_start and i < len(modified):
        if modified[i] == "(":
            j = modified.find(")", i)
            i = j + 1 if j != -1 else i + 1
            continue
        aa_idx += 1
        i += 1
    return aa_idx if aa_idx > 0 else None
