"""Supplementary tables for the final report."""

from __future__ import annotations

import json

import polars as pl

from megatensor.paths import EXPORTS, FIGURES, ROOT
from megatensor.store import PRIDE_STORE


def write_supplementary() -> dict[str, str]:
    EXPORTS.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    picks = FIGURES / "pride_glyco_picks.csv"
    if picks.is_file():
        df = pl.read_csv(picks).unique(subset=["accession"]).select(
            "accession",
            "title",
            "countries",
            "instrument_names",
            "software_names",
            "engine_kind",
            "likely_conditioned",
            "pick_note",
        )
        p = EXPORTS / "table_s1_pride_deposits.csv"
        df.write_csv(p)
        paths["table_s1"] = str(p)

    if PRIDE_STORE.summary_path.is_file():
        summary = json.loads(PRIDE_STORE.summary_path.read_text())
        pxds = summary.get("pxds", summary.get("datasets", []))
        spread = FIGURES / "pride_engine_spread.csv"
        if spread.is_file():
            eng = pl.read_csv(spread).filter(pl.col("dataset_id").is_in(pxds))
            p = EXPORTS / "table_s1_pride_stats.csv"
            eng.write_csv(p)
            paths["table_s1_stats"] = str(p)

    return paths
