"""§4 observation schema contract — column names and Polars dtypes."""

from __future__ import annotations

import polars as pl

# Observation columns emitted by every adapter (long-form, one row per metric).
OBSERVATION_COLUMNS: list[str] = [
    "dataset_id",
    "source_engine",
    "source_file",
    "file_checksum",
    "protein_id_raw",
    "protein_id_type",
    "isoform_raw",
    "residue_pos_raw",
    "residue_aa",
    "ptm_label",
    "ptm_unimod",
    "loc_score",
    "loc_method",
    "loc_is_ambiguous",
    "cond_cell_line",
    "cond_tissue",
    "cond_treatment",
    "cond_timepoint",
    "cond_fraction",
    "cond_replicate",
    "cond_replicate_type",
    "acq_ms_mode",
    "acq_msn_level",
    "acq_gradient_min",
    "acq_collision",
    "acq_enrichment",
    "inst_vendor",
    "inst_model",
    "inst_ms_cv",
    "prov_pxd",
    "prov_lab_pi",
    "prov_country",
    "prov_doi",
    "prov_search_software",
    "prov_search_version",
    "metric_name",
    "metric_value",
    "metric_norm_state",
    "metric_unit",
    "qc_flags",
]

# Registry-resolved columns added after §7.
RESOLVED_COLUMNS: list[str] = [
    "protein_acc",
    "isoform",
    "residue_pos",
    "identity_id",
    "ptm_id",
    "condition_id",
    "acquisition_id",
    "instrument_id",
    "provenance_id",
    "set_uid",
]


def empty_observations() -> pl.DataFrame:
    return pl.DataFrame({c: [] for c in OBSERVATION_COLUMNS}).cast(
        {
            "residue_pos_raw": pl.Int64,
            "loc_score": pl.Float64,
            "loc_is_ambiguous": pl.Boolean,
            "acq_gradient_min": pl.Float64,
            "metric_value": pl.Float64,
        }
    )
