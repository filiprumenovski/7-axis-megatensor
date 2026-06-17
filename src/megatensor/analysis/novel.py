"""Novel vs canon-shared site characterization."""

from __future__ import annotations

import polars as pl

from megatensor.analysis.common import human_site_keys, load_identity, load_site_index
from megatensor.store import UNION_STORE
from megatensor.viz.completeness import is_human_uniprot


def site_gene_lists() -> dict[str, list[str]]:
    idx = load_site_index()
    enrich_path = UNION_STORE.enrichment / "site_features.parquet"
    if not enrich_path.is_file():
        return {}

    enrich = pl.read_parquet(enrich_path)
    enrich = enrich.with_columns(
        pl.concat_str(
            [pl.col("protein_acc"), pl.col("residue_pos").cast(pl.Utf8), pl.col("residue_aa")],
            separator=":",
        ).alias("site_key")
    )

    pride_only_keys = set(
        idx.filter(pl.col("n_layers") == 1, pl.col("layers").list.contains("pride"))["site_key"].to_list()
    )
    shared_keys = set(idx.filter(pl.col("n_layers") >= 2)["site_key"].to_list())

    pride_only_keys = set(human_site_keys(list(pride_only_keys)))
    shared_keys = set(human_site_keys(list(shared_keys)))

    def genes(keys: set[str]) -> list[str]:
        sub = enrich.filter(pl.col("site_key").is_in(list(keys)))
        return [g for g in sub["gene_symbol"].drop_nulls().unique().sort().to_list() if g]

    return {
        "pride_novel": genes(pride_only_keys),
        "canon_shared": genes(shared_keys),
    }


def novel_summary() -> pl.DataFrame:
    idx = load_site_index()
    pride_id = load_identity("pride").filter(~pl.col("protein_level_only"))
    keys = idx.with_columns(
        pl.col("site_key").str.split(":").list.get(0).alias("protein_acc"),
    )
    human = keys.filter(pl.col("protein_acc").map_elements(lambda a: is_human_uniprot(a), return_dtype=pl.Boolean))

    return pl.DataFrame(
        {
            "class": ["canon_only", "shared", "pride_novel_human"],
            "n_sites": [
                human.filter(pl.col("n_layers") == 1, pl.col("layers").list.contains("canon")).height,
                human.filter(pl.col("n_layers") >= 2).height,
                human.filter(pl.col("n_layers") == 1, pl.col("layers").list.contains("pride")).height,
            ],
        }
    )
