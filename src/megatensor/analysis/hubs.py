"""Protein-level replication hubs with gene symbols."""

from __future__ import annotations

import polars as pl

from megatensor.enrich.uniprot import fetch_batch, gene_symbol


def hubs_with_genes(hubs: pl.DataFrame) -> pl.DataFrame:
    if hubs.is_empty():
        return hubs
    accs = hubs["protein_acc"].unique().to_list()
    entries = fetch_batch(accs)
    symbols = {acc: gene_symbol(entries.get(acc, {})) for acc in accs}
    return hubs.with_columns(
        pl.col("protein_acc")
        .map_elements(lambda a: symbols.get(a), return_dtype=pl.Utf8)
        .alias("gene_symbol")
    )
