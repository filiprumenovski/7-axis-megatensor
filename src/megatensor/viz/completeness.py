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
    "provenance": ["dataset_id", "prov_pxd", "source_engine"],
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
    acc_col = "protein_acc" if "protein_acc" in df.columns else "protein_id_raw"
    pos_col = "residue_pos" if "residue_pos" in df.columns else "residue_pos_raw"
    sites = df.select(
        "dataset_id",
        pl.concat_str(
            [pl.col(acc_col), pl.col(pos_col).cast(pl.Utf8), pl.col("residue_aa")],
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


def pride_engine_spread(df: pl.DataFrame) -> pl.DataFrame:
    """Per-PXD engine / instrument / country coverage for the PRIDE narrative."""
    spread = (
        df.group_by("dataset_id")
        .agg(
            pl.col("source_engine").drop_nulls().unique().sort().alias("engines"),
            pl.col("inst_model").drop_nulls().unique().sort().alias("instruments"),
            pl.col("prov_country").drop_nulls().unique().sort().alias("countries"),
            pl.col("acq_ms_mode").drop_nulls().unique().sort().alias("ms_modes"),
            pl.struct(["protein_id_raw", "residue_pos_raw", "residue_aa"]).n_unique().alias("unique_sites"),
            pl.len().alias("observation_rows"),
        )
        .sort("dataset_id")
    )
    return spread.with_columns(
        pl.col("engines").list.join("|"),
        pl.col("instruments").list.join("|"),
        pl.col("countries").list.join("|"),
        pl.col("ms_modes").list.join("|"),
    )


def site_index(sites: pl.DataFrame) -> pl.DataFrame:
    return (
        sites.filter(
            pl.col("site_key").is_not_null(),
            pl.col("site_key").str.contains(r"^[^:]+:\d+:[ST]$"),
        )
        .group_by("site_key")
        .agg(
            pl.col("layer").unique().sort().alias("layers"),
            pl.col("dataset_id").unique().sort().alias("datasets"),
            pl.len().alias("set_hits"),
        )
        .with_columns(pl.col("layers").list.len().alias("n_layers"))
        .sort("n_layers", descending=True)
    )


def cross_layer_overlap(canon: pl.DataFrame, pride: pl.DataFrame) -> pl.DataFrame:
    canon_keys = set(canon["site_key"].unique().to_list())
    pride_keys = set(pride["site_key"].unique().to_list())
    shared = canon_keys & pride_keys
    rows = [
        {"metric": "canon_sites", "value": len(canon_keys)},
        {"metric": "pride_sites", "value": len(pride_keys)},
        {"metric": "shared_sites", "value": len(shared)},
        {"metric": "canon_only", "value": len(canon_keys - pride_keys)},
        {"metric": "pride_only", "value": len(pride_keys - canon_keys)},
        {
            "metric": "jaccard",
            "value": round(len(shared) / len(canon_keys | pride_keys), 4) if canon_keys | pride_keys else 0.0,
        },
    ]
    return pl.DataFrame(rows)


def is_human_uniprot(acc: str | None) -> bool:
    if not acc:
        return False
    if acc.startswith("Os") or acc.startswith("ENST"):
        return False
    return acc[0].isalpha() and acc.isalnum() and 6 <= len(acc) <= 10


def human_pride_sites(pride_sites: pl.DataFrame) -> pl.DataFrame:
    """Exclude plant locus IDs mis-tagged as UniProt (e.g. rice PXD036527)."""
    return pride_sites.filter(
        pl.col("protein_acc").map_elements(is_human_uniprot, return_dtype=pl.Boolean)
    )
