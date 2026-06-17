"""PRIDE bronze metadata for provenance axis enrichment."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import polars as pl

from megatensor.paths import PRIDE
from megatensor.pipeline.pride_discover import find_snapshot_date


@dataclass(frozen=True)
class PrideProjectMeta:
    pxd: str
    title: str
    country: str | None
    instruments: str | None
    software: str | None
    inst_vendor: str | None
    inst_model: str | None
    doi: str | None = None
    lab_pi: str | None = None


def _split_instrument(name: str | None) -> tuple[str | None, str | None]:
    if not name:
        return None, None
    vendors = ("Thermo", "Bruker", "Waters", "AB Sciex", "Agilent", "Shimadzu")
    for v in vendors:
        if name.startswith(v):
            return v, name
    if "Orbitrap" in name or "LTQ" in name or "Q Exactive" in name or "Exploris" in name or "Astral" in name:
        return "Thermo", name
    if "timsTOF" in name:
        return "Bruker", name
    return None, name


@lru_cache(maxsize=1)
def load_project_meta(pride_root: str | None = None) -> dict[str, PrideProjectMeta]:
    root = Path(pride_root) if pride_root else PRIDE
    snap = find_snapshot_date(root)
    if not snap:
        return {}
    import duckdb

    glob = (root / "bronze" / "projects" / f"snapshot_date={snap}" / "*.parquet").as_posix()
    rows = duckdb.connect().execute(
        f"""
        SELECT
          accession,
          title,
          countries[1] AS country,
          list_transform(instruments, x -> x.name)[1] AS instrument,
          array_to_string(list_transform(softwares, x -> x.name), '|') AS software,
          doi,
          list_transform(labPIs, x -> x.name)[1] AS lab_pi
        FROM read_parquet('{glob}')
        """
    ).fetchall()
    out: dict[str, PrideProjectMeta] = {}
    for pxd, title, country, inst, sw, doi, lab_pi in rows:
        vendor, model = _split_instrument(inst)
        out[pxd] = PrideProjectMeta(
            pxd=pxd,
            title=title or "",
            country=country,
            instruments=inst,
            software=sw,
            inst_vendor=vendor,
            inst_model=model,
            doi=doi,
            lab_pi=lab_pi,
        )
    return out


def meta_for(pxd: str, pride_root: str | None = None) -> PrideProjectMeta:
    return load_project_meta(pride_root).get(
        pxd,
        PrideProjectMeta(pxd, "", None, None, None, None, None, None, None),
    )
