# 7-Axis Megatensor: Technical Report

## Summary Table

| Metric | Canonical Databases | PRIDE Experimental Deposits |
|---|---:|---:|
| Observation rows | 85,833 | 91,376 |
| Unique sites | 43,853 | 13,567 |
| SET records | 81,486 | 36,737 |

- **Cross-layer recovery:** 4,376 shared sites; 9,191 PRIDE-only; 36,818 canon-only.
- **Human PRIDE subset** (rice PXD036527 excluded): 4,376 shared; 8,925 PRIDE-only novel sites.
- **Deposit diversity:** Twelve PRIDE projects spanning MaxQuant, Proteome Discoverer, and mzTab across US and Chinese laboratories on Orbitrap Elite, Fusion, and Lumos instruments.

---

## Abstract

We developed a sparse, append-only Site Event Tensor (SET) megatensor for O-GlcNAc proteomics. The system harmonizes two canonical reference libraries (O-GlcNAc Database, O-GlcNAcAtlas 4.0) and twelve PRIDE experimental deposits into a shared seven-axis contract. This structure enables cross-layer site queries, provenance slicing, and machine learning exports without pairwise manual table reconciliation.

---

## Seven-Axis SET Schema

- **Layer A (Site Identity):** Protein accession, residue position, modified amino acid (Ser or Thr), and modification class (O-GlcNAc / UniMod:43).
- **Layer B (Experimental Context):** Quantification, biological condition, acquisition mode, instrument model, and repository provenance. Localization confidence scores are stored as payloads on the SET rather than as identity coordinates.

---

## Analysis Figures

- **Figure 0:** `figures/analysis_megatensor_impact.pdf` (Impact summary: flat file burden vs. megatensor ingest).
- **Figure 1:** `figures/analysis_replication_tiers.pdf` (Replication tiers: PRIDE-only, multi-PXD, canon-intersected, triangulated).
- **Figure 2:** `figures/analysis_protein_hubs.pdf` (Top protein hubs: HCFC1, NUP214, NUP98, OGT).
- **Figure 3 & 3b:** `figures/analysis_silac_ma.pdf` and `figures/analysis_silac_triangulation.pdf` (SILAC M-A plot and triangulation overlay).
- **Figure 4, 4b, 4c & 7:** Concordance suite (`figures/analysis_concordance_scatter.pdf`, `figures/analysis_concordance_context.pdf`, `figures/analysis_ogt_concordance.pdf`, `figures/analysis_concordance_heatmap.pdf`).
- **Figure 5:** `figures/analysis_chemoproteomics.pdf` (Chemoproteomic probe selectivity).
- **Figure 6 & 6b:** `figures/analysis_bap1ko_tissue.pdf` and `figures/analysis_bap1ko_brain_liver.pdf` (BAP1KO tissue burden and brain vs. liver contrast).
- **Figure 8:** `figures/analysis_gsea_contrast.pdf` (GO Biological Process pathway enrichment contrast).
- **Figure 9:** `figures/analysis_evidence_ladder.pdf` (Composite evidence score ranking).
- **Figure 10:** `figures/analysis_triangulated_heatmap.pdf` (Triangulated site observation intensity matrix).

---

## SQL Query Interface

```sql
-- Query observations across studies for a specific site
SELECT dataset_id, cond_treatment, metric_name, metric_value, inst_model, prov_country
FROM read_parquet('megatensor/pride/staging/observations.parquet')
WHERE protein_id_raw = 'Q6ZU65' AND residue_pos_raw = 1003;
```

Full query suite available in `queries/queries.sql`.

---

## Reproducibility

```bash
just setup && just download && just canon
just pride-discover && just pride-download && just pride-tensorize
just union && just analyze && just enrich && just export
duckdb < queries/queries.sql
```

---

## References

1. Wulff-Fuentes E, et al. The human O-GlcNAcome database and meta-analysis. *Scientific Data*. 2021;8:25.
2. Ma J, et al. O-GlcNAcAtlas: a database of experimentally identified O-GlcNAc sites and proteins. *Glycobiology*. 2021;31(7):719-723.
3. Hou C, et al. O-GlcNAcAtlas 4.0: An updated protein O-GlcNAcylation database with site-specific quantification. *Journal of Molecular Biology*. 2025;437(15):169033.
4. PRIDE Archive. European Bioinformatics Institute. <https://www.ebi.ac.uk/pride/>.
5. Rumenovski F. *pride-ingest*: Reproducible PRIDE metadata ingestion. <https://github.com/filiprumenovski/pride-ingest>.
