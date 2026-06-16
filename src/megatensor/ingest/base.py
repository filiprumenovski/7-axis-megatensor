"""Adapter interface and shared helpers."""

from __future__ import annotations

import hashlib
from abc import ABC, abstractmethod
from pathlib import Path

import polars as pl

from megatensor.hash_utils import PTM_UNIMOD
from megatensor.schema import OBSERVATION_COLUMNS


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_site_token(token: str) -> tuple[str | None, int | None]:
    """Parse 'S493' or 'T12' -> (aa, pos)."""
    token = token.strip().upper()
    if len(token) < 2:
        return None, None
    aa, pos_s = token[0], token[1:]
    if aa not in ("S", "T") or not pos_s.isdigit():
        return None, None
    return aa, int(pos_s)


def melt_sites_column(
    df: pl.DataFrame,
    sites_col: str,
    *,
    dataset_id: str,
    source_engine: str,
    source_file: str,
    file_checksum: str,
    protein_col: str,
    protein_id_type: str,
    loc_method: str,
    loc_is_ambiguous: bool,
    prov_doi: str | None = None,
) -> pl.DataFrame:
    """Explode semicolon-separated site tokens (MCW style) into observation rows."""
    base = df.select(
        pl.col(protein_col).alias("protein_id_raw"),
        pl.lit(protein_id_type).alias("protein_id_type"),
        pl.lit(None).cast(pl.Utf8).alias("isoform_raw"),
        pl.col(sites_col).alias("_sites"),
    ).filter(pl.col("_sites").is_not_null() & (pl.col("_sites") != ""))

    rows: list[pl.DataFrame] = []
    for row in base.iter_rows(named=True):
        sites = str(row["_sites"]).split(";")
        for site in sites:
            aa, pos = parse_site_token(site)
            if aa is None:
                continue
            rows.append(
                pl.DataFrame(
                    {
                        "dataset_id": [dataset_id],
                        "source_engine": [source_engine],
                        "source_file": [source_file],
                        "file_checksum": [file_checksum],
                        "protein_id_raw": [row["protein_id_raw"]],
                        "protein_id_type": [protein_id_type],
                        "isoform_raw": [None],
                        "residue_pos_raw": [pos],
                        "residue_aa": [aa],
                        "ptm_label": ["O-GlcNAc"],
                        "ptm_unimod": [PTM_UNIMOD],
                        "loc_score": [None],
                        "loc_method": [loc_method],
                        "loc_is_ambiguous": [loc_is_ambiguous],
                        "cond_cell_line": [None],
                        "cond_tissue": [None],
                        "cond_treatment": [None],
                        "cond_timepoint": [None],
                        "cond_fraction": [None],
                        "cond_replicate": [None],
                        "cond_replicate_type": [None],
                        "acq_ms_mode": [None],
                        "acq_msn_level": [None],
                        "acq_gradient_min": [None],
                        "acq_collision": [None],
                        "acq_enrichment": [None],
                        "inst_vendor": [None],
                        "inst_model": [None],
                        "inst_ms_cv": [None],
                        "prov_pxd": [None],
                        "prov_lab_pi": [None],
                        "prov_country": [None],
                        "prov_doi": [prov_doi],
                        "prov_search_software": [None],
                        "prov_search_version": [None],
                        "metric_name": ["spectral_count"],
                        "metric_value": [1.0],
                        "metric_norm_state": ["curated"],
                        "metric_unit": [None],
                        "qc_flags": [None],
                    }
                )
            )
    if not rows:
        return pl.DataFrame({c: [] for c in OBSERVATION_COLUMNS})
    return pl.concat(rows).select(OBSERVATION_COLUMNS)


class CanonAdapter(ABC):
    dataset_id: str
    source_engine: str

    @abstractmethod
    def parse(self, path: Path) -> pl.DataFrame:
        """Return long-form observation rows."""
