# 7-Axis Megatensor

**Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model**

UROP project (Filip Rumenovski / Dr. Charlie Fehl, WSU). Sparse, append-only Site Event Tensors (SETs) unioned into a Megatensor queryable via DuckDB.

Build doctrine: [`7_AXIS_MEGATENSOR_BUILD.md`](7_AXIS_MEGATENSOR_BUILD.md).
Figure style: [`FIGURES.md`](FIGURES.md). UROP final report: [`biorxiv.md`](biorxiv.md).

## Quick start

```bash
just setup      # venv + pip install -e .
just download   # bulk CSV: O-GlcNAc DB (MCW) + O-GlcNAcAtlas Dataset-I/II
just canon      # Phase 0: canon adapters -> observations -> SETs
pip install -e ".[pride]"   # optional: live re-ingest only
just unpack-pride           # once: unpack pride_snapshot_parquets_*.tar.gz
just pride-discover         # query local snapshot -> glyco candidates + ranked picks
just pride-download         # Aspera pull of curated result tables
just pride-tensorize        # Phase 3: PRIDE -> isolated experimental tensor
just union && just analyze && just analysis-figures
just enrich && just export && just biorxiv   # UROP report bundle
just publish                               # all of the above
duckdb < queries/queries.sql
```

## Data layout

```text
data/canon/           # bulk reference CSVs (gitignored)
data/pride/           # pride-ingest bronze/silver (gitignored)
megatensor/           # on-disk SET store (gitignored)
  registry/           # identity, ptm, condition, ... dim tables
  sets/set_coordinates/
  metrics/set_metrics/
views/megatensor.sql  # DuckDB union view
```

## Citations

- Wulff-Fuentes et al. 2021, PMID [33479245](https://pubmed.ncbi.nlm.nih.gov/33479245) (O-GlcNAc Database, MCW)
- Ma et al. 2021, Glycobiology; Hou et al. 2025, J Mol Biol (O-GlcNAcAtlas)
- PRIDE Archive via [`pride-ingest`](https://github.com/filiprumenovski/pride-ingest)

## Scope

This run produces ML-*ready* exports only. No model training.
