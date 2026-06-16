"""Figure A: per-source axis completeness."""

from __future__ import annotations

import polars as pl

AXIS_COLS = {
    "identity": ["protein_acc", "residue_pos", "residue_aa"],
    "ptm": ["ptm_unimod"],
    "quant": ["metric_value"],
    "condition": ["cond_cell_line", "cond_tissue", "cond_treatment", "cond_timepoint"],
    "acquisition": ["acq_ms_mode", "acq_msn_level", "acq_enrichment"],
    "instrument": ["inst_vendor", "inst_model"],
    "provenance": ["dataset_id", "prov_doi", "source_engine"],
}


def _axis_filled(df: pl.DataFrame, cols: list[str]) -> float:
    if not cols:
        return 0.0
    present = df.select([pl.col(c).is_not_null().alias(c) for c in cols if c in df.columns])
    if present.width == 0:
        return 0.0
    row_ok = present.select(pl.all_horizontal(pl.all()).alias("ok"))["ok"]
    return float(row_ok.mean()) if row_ok.len() else 0.0


def axis_completeness(df: pl.DataFrame) -> pl.DataFrame:
    rows = []
    for dataset_id in df["dataset_id"].unique().sort().to_list():
        part = df.filter(pl.col("dataset_id") == dataset_id)
        row = {
            "dataset_id": dataset_id,
            "set_count": part["set_uid"].n_unique() if "set_uid" in part.columns else part.height,
            "protein_count": part["protein_acc"].n_unique() if "protein_acc" in part.columns else 0,
            "site_count": part.select(["protein_acc", "residue_pos", "residue_aa"]).unique().height,
            "metric_count": part.height,
        }
        for axis, cols in AXIS_COLS.items():
            row[f"pct_{axis}"] = round(100 * _axis_filled(part, cols), 2)
        rows.append(row)
    return pl.DataFrame(rows)


def canon_overlap(df: pl.DataFrame) -> pl.DataFrame:
    """Site keys per canon source for UpSet / overlap stats."""
    sites = df.select(
        "dataset_id",
        pl.concat_str(
            [pl.col("protein_acc"), pl.col("residue_pos").cast(pl.Utf8), pl.col("residue_aa")],
            separator=":",
        ).alias("site_key"),
    ).unique()
    sources = sites["dataset_id"].unique().sort().to_list()
    by_source = {s: set(sites.filter(pl.col("dataset_id") == s)["site_key"].to_list()) for s in sources}
    rows = []
    for a in sources:
        for b in sources:
            if a >= b:
                continue
            inter = len(by_source[a] & by_source[b])
            rows.append(
                {
                    "source_a": a,
                    "source_b": b,
                    "overlap_sites": inter,
                    "only_a": len(by_source[a] - by_source[b]),
                    "only_b": len(by_source[b] - by_source[a]),
                }
            )
    return pl.DataFrame(rows)
