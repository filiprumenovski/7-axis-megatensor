# 7-Axis Megatensor

**Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model**

Undergraduate Research Opportunities Program (UROP) project by Filip Rumenovski (Faculty Mentor: Dr. Charlie Fehl, Department of Chemistry, Wayne State University).

Sparse, append-only Site Event Tensors (SETs) unified into a queryable Megatensor via DuckDB and Polars.

- Technical build specification: [`7_AXIS_MEGATENSOR_BUILD.md`](7_AXIS_MEGATENSOR_BUILD.md)
- Figure style standards: [`FIGURES.md`](FIGURES.md)
- Architectural rationale: [`WHY_MEGATENSOR.md`](WHY_MEGATENSOR.md)
- UROP final report / preprint draft: [`biorxiv.md`](biorxiv.md)
- Technical summary: [`report.md`](report.md)

## Quick Start

```bash
just setup      # virtual environment and dependencies
just download   # fetch reference databases: O-GlcNAc DB (MCW) + O-GlcNAcAtlas
just canon      # build canonical reference tensors
just pride-discover   # scan PRIDE snapshot for glycoproteomics datasets
just pride-download   # download curated PRIDE result tables
just pride-tensorize  # tensorize PRIDE experimental tables
just union      # merge canonical and PRIDE layers
just analyze    # calculate replication tiers, SILAC ratios, and concordance
just enrich     # attach sequence window and pathway annotations
just export     # export ML-ready feature matrices
duckdb < queries/queries.sql
```

## Data Architecture

```text
data/canon/           # reference database CSVs (gitignored)
data/pride/           # PRIDE result tables (gitignored)
megatensor/           # on-disk SET Parquet store (gitignored)
  registry/           # protein identity, residue coordinates
  sets/set_coordinates/
  metrics/set_metrics/
views/megatensor.sql  # DuckDB relational union view
exports/              # ML-ready site-by-condition and site-by-feature tensors
```

## Scope and Boundaries

The pipeline generates structured, standardized Parquet and NumPy matrices for machine learning and comparative analysis. Model training is out of scope for this pipeline.

## Citations

- Wulff-Fuentes et al. 2021, *Scientific Data*, PMID [33479245](https://pubmed.ncbi.nlm.nih.gov/33479245) (O-GlcNAc Database, MCW)
- Ma et al. 2021, *Glycobiology*; Hou et al. 2025, *J Mol Biol* (O-GlcNAcAtlas)
- PRIDE Archive via [`pride-ingest`](https://github.com/filiprumenovski/pride-ingest)
