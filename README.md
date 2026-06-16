# FEHL Megatensor

**Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model**

UROP project (Filip Rumenovski / Dr. Charlie Fehl, WSU). Sparse, append-only Site Event Tensors (SETs) unioned into a Megatensor queryable via DuckDB.

Build doctrine: [`FEHL_MEGATENSOR_BUILD.md`](FEHL_MEGATENSOR_BUILD.md).

## Quick start

```bash
just setup      # venv + pip install -e .
just download   # bulk CSV: O-GlcNAc DB (MCW) + O-GlcNAcAtlas Dataset-I/II
just canon      # Phase 0: canon adapters -> observations -> SETs
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
