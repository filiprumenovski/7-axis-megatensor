# Why the Megatensor Matters

*Technical justification and quantitative evidence from pipeline outputs (`figures/megatensor_importance.json`).*

---

## 1. The Core Problem

O-GlcNAc proteomics datasets are deposited in conflicting tabular formats. Search engines and laboratories mix biological site identity (what was modified) with experimental context (instrument settings, chemical probes, tissue types, and search software). As a result, cross-study replication and machine learning require extensive manual table formatting.

---

## 2. The Architectural Solution

We designed a sparse **Site Event Tensor (SET)** representation:
- Seven orthogonal axes separate biological identity from experimental context.
- Canonical reference databases and experimental PRIDE deposits are unified into an append-only store.
- Data queries are executed directly through DuckDB and Polars.
- Standardized matrices are exported for machine learning.

Ingestion uses thin per-format adapters. Downstream harmonization into the shared coordinate system is automated.

---

## 3. Quantitative Payoffs

### 3.1 Structural Interoperability

| Metric | Flat-File Approach | Megatensor Architecture |
|---|---|---|
| Pairwise conversions (12 PRIDE studies) | 66 manual conversion scripts | 12 one-time adapters |
| Site alignment | Manual string matching per comparison | Single join on `site_key` |
| Reference vs. repository overlap | Separate ad hoc scripts | 4,376 shared sites automatically |
| Replicated sites | Unstructured inspection | 353 triangulated sites (canon + >=2 PXDs) |

### 3.2 Context Preservation

A single stable identifier (`site_key`) connects observations across experimental conditions:
- **393** SILAC-quantified sites (heavy/light condition axis)
- **811** brain-liver site pairs (tissue condition axis)
- **4,622** chemoproteomic sites (probe chemistry axis)
- **82** sites that are both SILAC-quantified and triangulated across studies

### 3.3 Evidence Ranking and Machine Learning

- Stratified replication tiers and hub connectivity rankings (HCFC1: 18 triangulated sites across 3 PXDs).
- Empirical concordance testing demonstrates that intra-study replicates correlate strongly (r = 0.953), while cross-laboratory uncalibrated measurements correlate weakly (r = 0.274).
- Exported tensors: `site_x_condition.parquet` and `site_x_features.parquet`.

---

## 4. Analytical Capabilities Enabled by the Model

1. Direct ranking of 353 triangulated sites across canonical databases and multiple PRIDE accessions without pairwise joins.
2. Slicing tissue differences (brain vs. liver across 811 sites) and SILAC ratios using the same primary keys.
3. Quantifying cross-laboratory and intra-laboratory intensity concordance across search engines and facilities.
4. Slicing all observations for specific sites across studies in single SQL queries (see `queries/queries.sql`).

---

## 5. Scope and Boundaries

- **Intensity calibration:** The model establishes identity-level interoperability. Cross-study intensity scales remain laboratory-specific and require future normalization.
- **Localization scoring:** Search engine scores (ptmRS, MaxQuant localization probability, mzTab scores) are preserved as raw values rather than converted to an arbitrary unified scale.
- **Chemical discrimination:** Modification definitions are anchored to UniMod:43. Isobaric discrimination is reserved for future chemistry axes.

---

## 6. Reproduction

```bash
just analyze && just analysis-figures && just export
```
