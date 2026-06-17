"""UniProt sequence + feature fetch with on-disk cache."""

from __future__ import annotations

import json
import re
import time
from pathlib import Path

import httpx
import structlog

from megatensor.paths import CACHE

log = structlog.get_logger()

UNIPROT_BATCH = "https://rest.uniprot.org/uniprotkb/accessions"
CACHE_DIR = CACHE / "uniprot"
CHUNK = 100
ACC_RE = re.compile(r"^[A-Z][A-Z0-9]{5,9}$")


def _normalize_acc(token: str | None) -> str | None:
    if not token:
        return None
    base = str(token).split("-")[0].split(".")[0].strip().upper()
    if base.startswith("UNKNOWN") or base.startswith("REV__"):
        return None
    return base if ACC_RE.match(base) else None


def _cache_path(accession: str) -> Path:
    return CACHE_DIR / f"{accession}.json"


def fetch_batch(accessions: list[str], *, pause_s: float = 0.2) -> dict[str, dict]:
    """Return {accession: uniprot_json} for valid UniProt accessions."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out: dict[str, dict] = {}
    pending: list[str] = []

    for acc in accessions:
        norm = _normalize_acc(acc)
        if not norm:
            continue
        cached = _cache_path(norm)
        if cached.is_file():
            try:
                out[norm] = json.loads(cached.read_text())
                continue
            except json.JSONDecodeError:
                cached.unlink(missing_ok=True)
        pending.append(norm)

    for i in range(0, len(pending), CHUNK):
        chunk = pending[i : i + CHUNK]
        try:
            resp = httpx.post(
                UNIPROT_BATCH,
                params={"accessions": ",".join(chunk), "format": "json"},
                timeout=120.0,
            )
            resp.raise_for_status()
            payload = resp.json()
        except Exception as exc:
            log.warning("uniprot_batch_failed", chunk=len(chunk), error=str(exc))
            time.sleep(pause_s)
            continue

        for entry in payload.get("results", []):
            acc = entry.get("primaryAccession")
            if not acc:
                continue
            out[acc] = entry
            _cache_path(acc).write_text(json.dumps(entry))
        time.sleep(pause_s)

    return out


def sequence_of(entry: dict) -> str:
    seq = entry.get("sequence", {})
    return str(seq.get("value") or "")


def gene_symbol(entry: dict) -> str | None:
    for g in entry.get("genes", []) or []:
        name = g.get("geneName", {}).get("value")
        if name:
            return name
    return None


def features_at_position(entry: dict, pos: int) -> list[dict]:
    hits: list[dict] = []
    for feat in entry.get("features", []) or []:
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None or end is None:
            continue
        if int(start) <= pos <= int(end):
            hits.append(
                {
                    "type": feat.get("type"),
                    "description": feat.get("description"),
                    "start": int(start),
                    "end": int(end),
                }
            )
    return hits
