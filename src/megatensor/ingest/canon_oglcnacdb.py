"""O-GlcNAc Database (MCW) bulk CSV adapter."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from megatensor.ingest.base import CanonAdapter, file_sha256, melt_sites_column
from megatensor.schema import OBSERVATION_COLUMNS


class OGlcnacDbAdapter(CanonAdapter):
    dataset_id = "oglcnacdb"
    source_engine = "oglcnacdb"

    def parse(self, path: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
        checksum = file_sha256(path)
        rel = str(path.name)
        df = pl.read_csv(path, infer_schema_length=5000, ignore_errors=True)

        # Localized sites: semicolon-separated in "oglcnac sites"
        sites_col = "oglcnac sites"
        if sites_col not in df.columns:
            raise ValueError(f"expected column {sites_col!r} in {path}")

        protein_col = "UniprotKB ID" if "UniprotKB ID" in df.columns else "Uniprot ID"
        obs = melt_sites_column(
            df,
            sites_col,
            dataset_id=self.dataset_id,
            source_engine=self.source_engine,
            source_file=rel,
            file_checksum=checksum,
            protein_col=protein_col,
            protein_id_type="uniprot_acc",
            loc_method="manual_curation",
            loc_is_ambiguous=False,
            prov_doi="10.1038/s41597-021-00810-5",
        )

        # Protein-level-only rows (no site) — tracked separately, no SET emitted
        protein_only = df.filter(
            pl.col(sites_col).is_null() | (pl.col(sites_col).str.strip_chars() == "")
        ).select(
            pl.col(protein_col).alias("protein_id_raw"),
            pl.lit("uniprot_acc").alias("protein_id_type"),
            pl.lit(True).alias("protein_level_only"),
        )
        return obs, protein_only  # type: ignore[return-value]
