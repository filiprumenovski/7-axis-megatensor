"""O-GlcNAcAtlas Dataset-I / Dataset-II bulk CSV adapters."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from megatensor.hash_utils import PTM_UNIMOD
from megatensor.ingest.base import CanonAdapter, file_sha256
from megatensor.schema import OBSERVATION_COLUMNS


def _atlas_rows(
    df: pl.DataFrame,
    *,
    dataset_id: str,
    source_file: str,
    file_checksum: str,
    loc_method: str,
    loc_is_ambiguous: bool,
) -> pl.DataFrame:
    base = df.with_columns(
        pl.lit(dataset_id).alias("dataset_id"),
        pl.lit("atlas").alias("source_engine"),
        pl.lit(source_file).alias("source_file"),
        pl.lit(file_checksum).alias("file_checksum"),
        pl.col("accession").alias("protein_id_raw"),
        pl.lit("uniprot_acc").alias("protein_id_type"),
        pl.lit(None).cast(pl.Utf8).alias("isoform_raw"),
        pl.col("position_in_protein").cast(pl.Int64).alias("residue_pos_raw"),
        pl.col("site_residue").str.to_uppercase().alias("residue_aa"),
        pl.lit("O-GlcNAc").alias("ptm_label"),
        pl.lit(PTM_UNIMOD).alias("ptm_unimod"),
        pl.lit(None).cast(pl.Float64).alias("loc_score"),
        pl.lit(loc_method).alias("loc_method"),
        pl.lit(loc_is_ambiguous).alias("loc_is_ambiguous"),
        pl.col("sample_type").alias("cond_cell_line"),
        pl.lit(None).cast(pl.Utf8).alias("cond_tissue"),
        pl.col("condition").alias("cond_treatment"),
        pl.lit(None).cast(pl.Utf8).alias("cond_timepoint"),
        pl.lit(None).cast(pl.Utf8).alias("cond_fraction"),
        pl.lit(None).cast(pl.Utf8).alias("cond_replicate"),
        pl.lit(None).cast(pl.Utf8).alias("cond_replicate_type"),
        pl.lit(None).cast(pl.Utf8).alias("acq_ms_mode"),
        pl.lit(None).cast(pl.Utf8).alias("acq_msn_level"),
        pl.lit(None).cast(pl.Float64).alias("acq_gradient_min"),
        pl.lit(None).cast(pl.Utf8).alias("acq_collision"),
        pl.lit(None).cast(pl.Utf8).alias("acq_enrichment"),
        pl.lit(None).cast(pl.Utf8).alias("inst_vendor"),
        pl.lit(None).cast(pl.Utf8).alias("inst_model"),
        pl.lit(None).cast(pl.Utf8).alias("inst_ms_cv"),
        pl.lit(None).cast(pl.Utf8).alias("prov_pxd"),
        pl.lit(None).cast(pl.Utf8).alias("prov_lab_pi"),
        pl.lit(None).cast(pl.Utf8).alias("prov_country"),
        pl.lit(None).cast(pl.Utf8).alias("prov_doi"),
        pl.lit(None).cast(pl.Utf8).alias("prov_search_software"),
        pl.lit(None).cast(pl.Utf8).alias("prov_search_version"),
    )

    # Curated quant events where log2Ratio present; else spectral_count=1
    has_quant = pl.col("log2Ratio").is_not_null() & (pl.col("log2Ratio") != "")
    quant = base.filter(has_quant).with_columns(
        pl.lit("quant_event").alias("metric_name"),
        pl.col("log2Ratio").cast(pl.Float64, strict=False).alias("metric_value"),
        pl.lit("log2").alias("metric_norm_state"),
        pl.lit(None).cast(pl.Utf8).alias("metric_unit"),
        pl.lit(None).cast(pl.Utf8).alias("qc_flags"),
    )
    site_only = base.filter(~has_quant).with_columns(
        pl.lit("spectral_count").alias("metric_name"),
        pl.lit(1.0).alias("metric_value"),
        pl.lit("curated").alias("metric_norm_state"),
        pl.lit(None).cast(pl.Utf8).alias("metric_unit"),
        pl.lit(None).cast(pl.Utf8).alias("qc_flags"),
    )
    return pl.concat([quant, site_only]).select(OBSERVATION_COLUMNS)


class AtlasDatasetIAdapter(CanonAdapter):
    dataset_id = "atlas_I"
    source_engine = "atlas"

    def parse(self, path: Path) -> pl.DataFrame:
        checksum = file_sha256(path)
        df = pl.read_csv(path, infer_schema_length=0, encoding="utf8-lossy").with_columns(
            pl.col("position_in_protein").str.extract(r"(\d+)", 1).cast(pl.Int64, strict=False),
            pl.col("position_in_peptide").cast(pl.Int64, strict=False),
        )
        return _atlas_rows(
            df,
            dataset_id=self.dataset_id,
            source_file=path.name,
            file_checksum=checksum,
            loc_method="atlas_unambiguous",
            loc_is_ambiguous=False,
        )


class AtlasDatasetIIAdapter(CanonAdapter):
    dataset_id = "atlas_II"
    source_engine = "atlas"

    def parse(self, path: Path) -> pl.DataFrame:
        checksum = file_sha256(path)
        df = pl.read_csv(path, infer_schema_length=0, encoding="utf8-lossy").with_columns(
            pl.col("position_in_protein").str.extract(r"(\d+)", 1).cast(pl.Int64, strict=False),
            pl.col("position_in_peptide").cast(pl.Int64, strict=False),
        )
        return _atlas_rows(
            df,
            dataset_id=self.dataset_id,
            source_file=path.name,
            file_checksum=checksum,
            loc_method="atlas_ambiguous",
            loc_is_ambiguous=True,
        )
