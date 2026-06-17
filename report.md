# FEHL Megatensor — Technical Report

## Panel headline

| | Canon | PRIDE |
|---|------:|------:|
| Observation rows | 85,833 | 91,376 |
| Unique sites | 43,853 | 13,567 |
| SETs | 81,486 | 36,737 |

**Cross-layer:** 4,376 shared sites · 9,191 PRIDE-only · 36,818 canon-only

**Human PRIDE** (rice PXD036527 excluded): 4,376 shared · 8,925 PRIDE-only novel vs canon

Twelve PRIDE deposits span **MaxQuant, Proteome Discoverer, and mzTab** across **US and China**, on **Orbitrap Elite / Fusion / Lumos** — all DDA.

## Abstract

We built a sparse, append-only **Site Event Tensor (SET)** megatensor for O-GlcNAc proteomics by harmonizing
two canonical reference libraries (O-GlcNAc Database + O-GlcNAcAtlas 4.0) and twelve PRIDE experimental
deposits into a shared 7-axis contract. The result supports cross-layer site queries, engine/instrument
provenance slicing, and ML-ready exports without pairwise manual reconciliation.

## Ontology (7 axes)

**Identity** (protein, position, S/T) and **PTM** (O-GlcNAc / UniMod:43) define the site entity.
**Quant**, **Condition**, **Acquisition**, **Instrument**, and **Provenance** encode experimental context.
Localization confidence (`loc_score`, `loc_method`) is payload on the SET, not part of the identity coordinate.

## Figures (PNG for slides)

- **figure_a:** `figures/figure_a_axis_completeness.pdf` (vector) · `figures/figure_a_axis_completeness.png` (preview)
- **figure_b:** `figures/figure_b_cross_layer_overlap.pdf` (vector) · `figures/figure_b_cross_layer_overlap.png` (preview)
- **figure_c:** `figures/figure_c_trajectory.pdf` (vector) · `figures/figure_c_trajectory.png` (preview)
- **pride_heterogeneity:** `figures/figure_pride_heterogeneity.pdf` (vector) · `figures/figure_pride_heterogeneity.png` (preview)
- **enrichment:** `figures/figure_enrichment_completeness.pdf` (vector) · `figures/figure_enrichment_completeness.png` (preview)

### Figure A — Axis completeness

Canon fills identity and reference provenance; PRIDE fills condition, acquisition, instrument, and engine-native quant metrics. CSV: `figures/axis_completeness_*.csv`.

### Figure B — Cross-layer interoperability

UpSet membership tables: `figures/upset_canon_membership.csv`, `figures/upset_pride_membership.csv`.
Human-filtered overlap: `figures/canon_vs_pride_overlap_human.json`.

### Figure C — Single-site trajectory

**Q6ZU65:1003:T — Light/Heavy SILAC in PXD039536**. Mean Light vs Heavy intensities in `figures/figure_c_trajectory.png`.
_SASA skipped (install `freesasa` or AlphaFold model unavailable)._

## Enrichment

| Feature | Coverage |
|---------|----------|
| Seq window | 73.66% |
| Domain/region | 82.36% |
| Disorder (metapredict) | 0.0% |
| Gene symbol | 85.24% |

_GSEA skipped (install `gseapy` or network unavailable)._

## Live demo query (panel)

```sql
-- One site, all PRIDE conditions that hit it
SELECT dataset_id, cond_treatment, metric_name, metric_value, inst_model, prov_country
FROM read_parquet('megatensor/pride/staging/observations.parquet')
WHERE protein_id_raw = 'Q6ZU65' AND residue_pos_raw = 1003;
```

Full query pack: `queries/queries.sql`

## Honest framing

Each source required a thin adapter mapping engine-native columns to the observation contract.
Harmonization **into the 7-axis space** is automatic downstream of that boundary; adapters are the
deliberate, auditable manual surface. Localization scores are **not** unified across engines (ptmRS vs MQ loc prob vs DIA-NN).

## Reproduction

```bash
just setup && just download && just canon
just pride-discover && just pride-download && just pride-tensorize
just union && just figures && just enrich && just export && just report
duckdb < queries/queries.sql
```

Install garnish: `pip install -e ".[enrich]"` (metapredict, gseapy, freesasa) and `pip install matplotlib`.

## Citations

- Wulff-Fuentes et al. 2021, *Sci Data* (O-GlcNAc Database)
- Ma et al. 2021; Hou et al. 2025 (O-GlcNAcAtlas)
- PRIDE Archive via pride-ingest
