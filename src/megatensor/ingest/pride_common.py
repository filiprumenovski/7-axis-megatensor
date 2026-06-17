"""Shared helpers for PRIDE file adapters."""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import polars as pl

from megatensor.hash_utils import PTM_UNIMOD
from megatensor.ingest.base import file_sha256
from megatensor.ingest.instruments import psi_ms_cv
from megatensor.ingest.pride_meta import PrideProjectMeta, meta_for
from megatensor.schema import OBSERVATION_COLUMNS

GLYCO_MOD_RE = re.compile(r"glcnac|hexnac|oglcna|o-glcnac|amtzhexnac|unimod:43|unimod:934|unimod:121", re.I)
UNIPROT_ACC_RE = re.compile(r"^[A-Z0-9]{6,10}")
SITE_MOD_RE = re.compile(r"([ST])(\d+)\(([^)]+)\)")


def read_table(path: Path) -> pl.DataFrame:
    if path.suffix.lower() == ".csv":
        return pl.read_csv(path, infer_schema_length=5000, ignore_errors=True, encoding="utf8-lossy")
    # PRIDE MaxQuant / PD exports are usually TSV; some PD tables are quoted TSV.
    try:
        return pl.read_csv(
            path,
            separator="\t",
            infer_schema_length=5000,
            ignore_errors=True,
            encoding="utf8-lossy",
            quote_char='"',
        )
    except Exception:
        return pl.read_csv(
            path,
            separator="\t",
            infer_schema_length=5000,
            ignore_errors=True,
            encoding="utf8-lossy",
        )


def parse_uniprot_acc(token: str | None) -> str | None:
    if not token:
        return None
    token = str(token).strip().split(";")[0].strip()
    if "|" in token:
        parts = token.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1].split("-")[0]
    base = token.split("-")[0].split(".")[0]
    if UNIPROT_ACC_RE.match(base):
        return base
    return base if base else None


def base_provenance(pxd: str, engine: str, path: Path, meta: PrideProjectMeta | None = None) -> dict:
    m = meta or meta_for(pxd)
    return {
        "dataset_id": pxd,
        "source_engine": engine,
        "source_file": f"{pxd}/{path.name}",
        "file_checksum": file_sha256(path),
        "prov_pxd": pxd,
        "prov_lab_pi": m.lab_pi,
        "prov_country": m.country,
        "prov_doi": m.doi,
        "prov_search_software": m.software,
        "prov_search_version": None,
        "inst_vendor": m.inst_vendor,
        "inst_model": m.inst_model,
        "inst_ms_cv": psi_ms_cv(m.inst_model),
        "ptm_label": "O-GlcNAc",
        "ptm_unimod": PTM_UNIMOD,
    }


def empty_obs() -> pl.DataFrame:
    return pl.DataFrame({c: [] for c in OBSERVATION_COLUMNS}).cast(
        {
            "residue_pos_raw": pl.Int64,
            "loc_score": pl.Float64,
            "loc_is_ambiguous": pl.Boolean,
            "acq_gradient_min": pl.Float64,
            "metric_value": pl.Float64,
        }
    )


def finalize_obs(df: pl.DataFrame) -> pl.DataFrame:
    if df.is_empty():
        return empty_obs()
    missing = [c for c in OBSERVATION_COLUMNS if c not in df.columns]
    for c in missing:
        df = df.with_columns(pl.lit(None).alias(c))
    return df.select(OBSERVATION_COLUMNS)


def condition_from_name(name: str) -> dict[str, str | None]:
    n = name.lower().replace("-", "_")
    out: dict[str, str | None] = {
        "cond_cell_line": None,
        "cond_tissue": None,
        "cond_treatment": None,
        "cond_timepoint": None,
        "cond_fraction": None,
        "cond_replicate": None,
        "cond_replicate_type": None,
    }
    if "293t" in n:
        out["cond_cell_line"] = "293T"
    if "liverbrain" in n or "liver_brain" in n:
        out["cond_tissue"] = "liver_brain"
    elif "brain" in n:
        out["cond_tissue"] = "brain"
    elif "liver" in n:
        out["cond_tissue"] = "liver"
    if "glycomics" in n:
        out["cond_fraction"] = "glycomics"
    elif "interactome" in n:
        out["cond_fraction"] = "interactome"
    elif "proteinexpression" in n or "protein_expression" in n:
        out["cond_fraction"] = "proteome"
    if "swissprot" in n:
        out["cond_replicate_type"] = "database_swissprot"
    elif "_full" in n or "full_" in n:
        out["cond_replicate_type"] = "database_full"
    if "light" in n:
        out["cond_treatment"] = "Light"
    if "heavy" in n:
        out["cond_treatment"] = "Heavy"
    if "insulin" in n:
        out["cond_treatment"] = "insulin"
    if "bap1ko" in n or "bap1_ko" in n:
        out["cond_treatment"] = "BAP1KO"
    elif "bap1" in n:
        out["cond_treatment"] = "BAP1KO"
    if re.search(r"\bwt\b", n):
        out["cond_treatment"] = "WT"
    if re.search(r"\bko\b", n) and "bap1" not in n:
        out["cond_treatment"] = "KO"
    for probe in ("pc", "dde", "dadps", "diazo"):
        if re.search(rf"\b{probe}\b", n) or f"_{probe}_" in n or n.startswith(f"{probe}_"):
            out["cond_treatment"] = probe.upper()
    return out


def collision_from_text(text: str | None) -> str | None:
    if not text:
        return None
    u = text.upper()
    if "ETHCD" in u or "ETD_HCD" in u or "HCD_ETD" in u:
        return "EThcD"
    if "HCD" in u:
        return "HCD"
    if "CID" in u:
        return "CID"
    if "ETD" in u:
        return "ETD"
    return None


def merge_conditions(*conds: dict[str, str | None]) -> dict[str, str | None]:
    merged: dict[str, str | None] = {
        "cond_cell_line": None,
        "cond_tissue": None,
        "cond_treatment": None,
        "cond_timepoint": None,
        "cond_fraction": None,
        "cond_replicate": None,
        "cond_replicate_type": None,
    }
    for cond in conds:
        for k, v in cond.items():
            if v and not merged.get(k):
                merged[k] = v
    return merged


def extract_from_zip(path: Path, suffix: str = ".txt") -> Path:
    """Extract first matching member beside the zip; return extracted path."""
    with zipfile.ZipFile(path) as zf:
        names = [n for n in zf.namelist() if n.lower().endswith(suffix.lower()) and not n.endswith("/")]
        preferred = [
            n
            for n in names
            if any(tok in n.lower() for tok in ("hexnac", "oglcnac", "o-glcnac", "og_site", "st_site"))
        ]
        member = preferred[0] if preferred else (names[0] if names else None)
        if not member:
            raise FileNotFoundError(f"no {suffix} in {path.name}")
        out = path.parent / Path(member).name
        if not out.exists() or out.stat().st_size == 0:
            out.write_bytes(zf.read(member))
        return out
