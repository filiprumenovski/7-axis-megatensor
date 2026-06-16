"""§8 controlled vocabularies and dimension tables."""

from __future__ import annotations

import polars as pl

from megatensor.hash_utils import hash_u64


def _dim_from_column(obs: pl.DataFrame, col: str, dim_name: str) -> tuple[pl.DataFrame, pl.Series]:
    """Build a dim table + id column from a nullable string axis column."""
    raw = obs.select(pl.col(col).fill_null("__NULL__").alias("token")).unique()
    dim = raw.with_columns(
        pl.col("token")
        .map_elements(lambda t: hash_u64(dim_name, t), return_dtype=pl.UInt64)
        .alias(f"{dim_name}_id"),
    )
    lookup = dict(zip(dim["token"].to_list(), dim[f"{dim_name}_id"].to_list(), strict=True))
    ids = obs[col].fill_null("__NULL__").map_elements(
        lambda t: lookup[t], return_dtype=pl.UInt64
    )
    return dim.rename({"token": f"{dim_name}_token"}), ids


def encode_context(obs: pl.DataFrame) -> tuple[pl.DataFrame, dict[str, pl.DataFrame]]:
    """Attach condition/acquisition/instrument/provenance ids; return dim tables."""
    dims: dict[str, pl.DataFrame] = {}
    out = obs

    condition_token = (
        pl.concat_str(
            [
                pl.col("cond_cell_line").fill_null(""),
                pl.col("cond_tissue").fill_null(""),
                pl.col("cond_treatment").fill_null(""),
                pl.col("cond_timepoint").fill_null(""),
                pl.col("cond_fraction").fill_null(""),
                pl.col("cond_replicate").fill_null(""),
            ],
            separator="|",
        )
        .alias("_cond")
    )
    out = out.with_columns(condition_token)
    cond_dim, cond_ids = _dim_from_column(out, "_cond", "condition")
    dims["condition_dim"] = cond_dim
    out = out.with_columns(cond_ids.alias("condition_id")).drop("_cond")

    acq_token = (
        pl.concat_str(
            [
                pl.col("acq_ms_mode").fill_null(""),
                pl.col("acq_msn_level").fill_null(""),
                pl.col("acq_collision").fill_null(""),
                pl.col("acq_enrichment").fill_null(""),
            ],
            separator="|",
        )
        .alias("_acq")
    )
    out = out.with_columns(acq_token)
    acq_dim, acq_ids = _dim_from_column(out, "_acq", "acquisition")
    dims["acquisition_dim"] = acq_dim
    out = out.with_columns(acq_ids.alias("acquisition_id")).drop("_acq")

    inst_token = (
        pl.concat_str(
            [pl.col("inst_vendor").fill_null(""), pl.col("inst_model").fill_null("")],
            separator="|",
        )
        .alias("_inst")
    )
    out = out.with_columns(inst_token)
    inst_dim, inst_ids = _dim_from_column(out, "_inst", "instrument")
    dims["instrument_dim"] = inst_dim
    out = out.with_columns(inst_ids.alias("instrument_id")).drop("_inst")

    prov_token = (
        pl.concat_str(
            [
                pl.col("dataset_id").fill_null(""),
                pl.col("prov_pxd").fill_null(""),
                pl.col("prov_country").fill_null(""),
                pl.col("prov_search_software").fill_null(""),
                pl.col("source_engine").fill_null(""),
            ],
            separator="|",
        )
        .alias("_prov")
    )
    out = out.with_columns(prov_token)
    prov_dim, prov_ids = _dim_from_column(out, "_prov", "provenance")
    dims["provenance_dim"] = prov_dim
    out = out.with_columns(prov_ids.alias("provenance_id")).drop("_prov")

    ptm_dim = pl.DataFrame({"ptm_unimod": ["UniMod:43"], "ptm_label": ["O-GlcNAc"], "ptm_id": [hash_u64("UniMod:43")]})
    dims["ptm_dim"] = ptm_dim

    return out, dims
