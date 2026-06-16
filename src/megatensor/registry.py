"""§7 identity normalization — UniProt accession + canonical coordinates."""

from __future__ import annotations

import polars as pl

from megatensor.hash_utils import PTM_ID, hash_u64


def resolve_identity(obs: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Resolve observations to registry fields.

  Canon bulk files already use UniProt accessions and protein positions.
  PRIDE adapters will plug in REST/ETL later.
    """
    resolved = obs.with_columns(
        pl.when(pl.col("protein_id_type").is_in(["uniprot_acc", "uniprot_id"]))
        .then(pl.col("protein_id_raw").str.split(".").list.first())
        .otherwise(pl.col("protein_id_raw"))
        .alias("protein_acc"),
        pl.col("isoform_raw").alias("isoform"),
        pl.col("residue_pos_raw").cast(pl.Int64).alias("residue_pos"),
    )

    # Flag non S/T and missing position
    resolved = resolved.with_columns(
        pl.when(~pl.col("residue_aa").is_in(["S", "T"]))
        .then(pl.concat_str([pl.col("qc_flags"), pl.lit("non_st")], separator="|"))
        .otherwise(pl.col("qc_flags"))
        .alias("qc_flags"),
    )

    identity_dim = (
        resolved.select(
            "protein_acc",
            "isoform",
            "residue_pos",
            "residue_aa",
        )
        .unique()
        .with_columns(
            pl.struct(["protein_acc", "isoform", "residue_pos", "residue_aa"])
            .map_elements(
                lambda r: hash_u64(r["protein_acc"], r["isoform"], r["residue_pos"], r["residue_aa"]),
                return_dtype=pl.UInt64,
            )
            .alias("identity_id"),
            pl.lit(False).alias("protein_level_only"),
        )
    )

    resolved = resolved.with_columns(
        pl.struct(["protein_acc", "isoform", "residue_pos", "residue_aa"])
        .map_elements(
            lambda r: hash_u64(r["protein_acc"], r["isoform"], r["residue_pos"], r["residue_aa"]),
            return_dtype=pl.UInt64,
        )
        .alias("identity_id"),
        pl.lit(PTM_ID).cast(pl.UInt64).alias("ptm_id"),
    )
    return resolved, identity_dim
