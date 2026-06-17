"""Cross-layer interoperability index (canon vs PRIDE site membership)."""

from __future__ import annotations

import json

import polars as pl
import structlog

from megatensor.paths import FIGURES, VIEWS
from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE
from megatensor.viz.completeness import cross_layer_overlap, human_pride_sites, site_index

log = structlog.get_logger()


def _load_layer_sites(store) -> pl.DataFrame:
    identity = pl.read_parquet(store.registry / "identity_dim.parquet")
    coords = pl.read_parquet(str(store.sets / "dataset_id=*/*.parquet"), hive_partitioning=True)
    joined = coords.join(identity, on="identity_id", how="left")
    return joined.with_columns(
        pl.lit(store.layer).alias("layer"),
        pl.concat_str(
            [
                pl.col("protein_acc"),
                pl.col("residue_pos").cast(pl.Utf8),
                pl.col("residue_aa"),
            ],
            separator=":",
        ).alias("site_key"),
    )


def run_union() -> dict:
    if not CANON_STORE.summary_path.exists():
        raise FileNotFoundError("canon store missing — run: just canon")
    if not PRIDE_STORE.summary_path.exists():
        raise FileNotFoundError("pride store missing — run: just pride-tensorize")

    canon_sites = _load_layer_sites(CANON_STORE)
    pride_sites = _load_layer_sites(PRIDE_STORE)
    all_sites = pl.concat([canon_sites, pride_sites], how="diagonal_relaxed")

    UNION_STORE.root.mkdir(parents=True, exist_ok=True)
    UNION_STORE.staging.mkdir(parents=True, exist_ok=True)

    index = site_index(all_sites)
    index.write_parquet(UNION_STORE.staging / "site_index.parquet")

    overlap = cross_layer_overlap(canon_sites, pride_sites)
    human_overlap = cross_layer_overlap(canon_sites, human_pride_sites(pride_sites))
    FIGURES.mkdir(parents=True, exist_ok=True)
    overlap.write_json(FIGURES / "canon_vs_pride_overlap.json")
    human_overlap.write_json(FIGURES / "canon_vs_pride_overlap_human.json")
    _write_union_views()

    summary = {
        "purpose": "cross_layer_interoperability_index",
        "canon_unique_sites": canon_sites.select("site_key").unique().height,
        "pride_unique_sites": pride_sites.select("site_key").unique().height,
        "shared_sites": int(overlap.filter(pl.col("metric") == "shared_sites")["value"][0])
        if overlap.height
        else 0,
        "canon_only_sites": int(overlap.filter(pl.col("metric") == "canon_only")["value"][0])
        if overlap.height
        else 0,
        "pride_only_sites": int(overlap.filter(pl.col("metric") == "pride_only")["value"][0])
        if overlap.height
        else 0,
        "human_pride_unique_sites": int(
            human_overlap.filter(pl.col("metric") == "pride_sites")["value"][0]
        )
        if human_overlap.height
        else 0,
        "human_shared_sites": int(human_overlap.filter(pl.col("metric") == "shared_sites")["value"][0])
        if human_overlap.height
        else 0,
        "human_pride_only_sites": int(
            human_overlap.filter(pl.col("metric") == "pride_only")["value"][0]
        )
        if human_overlap.height
        else 0,
    }
    UNION_STORE.summary_path.write_text(json.dumps(summary, indent=2))
    log.info("union_index_complete", **summary)
    print(json.dumps(summary, indent=2))
    return summary


def _write_union_views() -> None:
    root = CANON_STORE.root.parent
    VIEWS.mkdir(parents=True, exist_ok=True)
    (VIEWS / "canon.sql").write_text(
        f"""-- Canon reference tensor (identity backbone demo)
SELECT c.*, i.protein_acc, i.residue_pos, i.residue_aa, m.metric_name, m.metric_value
FROM read_parquet('{root}/canon/sets/set_coordinates/dataset_id=*/*.parquet', hive_partitioning := true) c
JOIN read_parquet('{root}/canon/registry/identity_dim.parquet') i USING (identity_id)
JOIN read_parquet('{root}/canon/metrics/set_metrics/dataset_id=*/*.parquet', hive_partitioning := true) m USING (set_uid, dataset_id);
"""
    )
    (VIEWS / "pride.sql").write_text(
        f"""-- PRIDE experimental tensor (rich context axes)
SELECT c.*, i.protein_acc, i.residue_pos, i.residue_aa,
       id.inst_model_token, pd.provenance_token, m.metric_name, m.metric_value
FROM read_parquet('{root}/pride/sets/set_coordinates/dataset_id=*/*.parquet', hive_partitioning := true) c
JOIN read_parquet('{root}/pride/registry/identity_dim.parquet') i USING (identity_id)
JOIN read_parquet('{root}/pride/registry/instrument_dim.parquet') id ON c.instrument_id = id.instrument_id
JOIN read_parquet('{root}/pride/registry/provenance_dim.parquet') pd ON c.provenance_id = pd.provenance_id
JOIN read_parquet('{root}/pride/metrics/set_metrics/dataset_id=*/*.parquet', hive_partitioning := true) m USING (set_uid, dataset_id);
"""
    )
    (VIEWS / "union_sites.sql").write_text(
        f"""-- Cross-layer site membership (Figure B interoperability)
SELECT * FROM read_parquet('{root}/union/staging/site_index.parquet');
"""
    )
