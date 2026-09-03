# 7-Axis Megatensor: Toward Interoperable O-GlcNAc Proteomics

**An open-source data architecture and harmonized compendium of 57,000+ protein O-GlcNAcylation sites across canonical databases and public PRIDE repositories.**

**Filip Rumenovski**¹ and **Charlie Fehl**¹*  
¹ *Department of Chemistry, College of Liberal Arts and Sciences, Wayne State University, Detroit, MI*  
* *Faculty Mentor & Principal Investigator: [Fehl Lab (Chemical Biology & Glycobiology)](https://fehllab.wayne.edu/)*  
* *Undergraduate Research Opportunities Program (UROP) Final Project*

---

## 1. Biological Motivation

Protein O-GlcNAcylation (the addition of a single $\beta$-N-acetylglucosamine sugar to serine and threonine residues) is a vital, nutrient-sensing post-translational modification. Unlike protein phosphorylation, which is governed by more than 500 distinct kinases, O-GlcNAcylation is orchestrated by a single pair of opposing enzymes:
- **OGT (O-GlcNAc transferase):** transfers GlcNAc from UDP-GlcNAc onto target proteins.
- **OGA (O-GlcNAcase):** hydrolyzes the glycosidic linkage to remove the modification.

Despite this centralized enzymatic control, O-GlcNAc modifies thousands of proteins across nuclear, cytosolic, and mitochondrial compartments, directly coordinating transcription, chromatin remodeling, nutrient sensing, and proteotoxic stress survival.

### The Data Fragmentation Problem
Over twenty years of mass spectrometry investigations have deposited thousands of O-GlcNAc datasets into public repositories like the PRIDE Archive and curated databases (such as the MCW O-GlcNAc Database and O-GlcNAcAtlas). However, these data remain largely non-interoperable:
- **Heterogeneous search outputs:** Different search algorithms (MaxQuant, Proteome Discoverer, MSFragger/FragPipe, DIA-NN) report conflicting column headers and incompatible localization probability formats.
- **Entangled metadata:** Biological site identity (`UniProt:position:AA`) and experimental variables (tissue type, chemical enrichment probe, SILAC channel, mass spectrometer model) are merged into flat spreadsheets.
- **Integration bottleneck:** Comparing twelve published datasets traditionally requires writing 66 separate pairwise conversion scripts, preventing cross-study validation and statistical meta-analysis.

The **7-Axis Megatensor** resolves this bottleneck by defining an append-only coordinate system that separates permanent biological site identity from variable experimental context.

---

## 2. Compendium Summary & Scale

The megatensor integrates two canonical reference databases and twelve diverse PRIDE experimental studies into a unified columnar store:

| Data Layer | Sources & Cohorts | Unique O-GlcNAc Sites | Total Observations | Search Software & Instrumentation |
|---|---|---:|---:|---|
| **Canonical Reference** | O-GlcNAc Database (MCW), O-GlcNAcAtlas 4.0 | 43,853 | 81,486 | Curated literature, multi-species |
| **PRIDE Public Archive** | 12 diverse PXDs (US & China cohorts) | 13,567 | 36,737 | MaxQuant, Proteome Discoverer, mzTab; Orbitrap Elite, Fusion, Lumos |
| **Harmonized Union** | **Integrated Compendium** | **53,044** | **118,223** | **7-Axis SET coordinate space** |

### Stratified Replication Tiers
- **4,376 sites** cross-intersect between raw PRIDE experimental tables and canonical reference databases without manual table formatting.
- **1,211 sites** are independently replicated in two or more PRIDE projects.
- **353 triangulated sites** represent highest-confidence biological benchmarks, supported by canonical literature and verified in at least two independent PRIDE deposits.

---

## 3. Key Scientific Findings

### 3.1 Regulatory Hubs Accumulate Replicated Sites
Replicated O-GlcNAc sites strongly concentrate on master regulatory proteins:
- **HCFC1 (Host Cell Factor 1, P51610):** The primary hub in the proteome, harboring **18 triangulated sites** verified across 3 independent PRIDE accessions.
- **Nuclear Pore Complex (NUP214, NUP98):** 14 and 12 triangulated sites, reflecting dense O-GlcNAc gating at the nuclear pore.
- **OGT Self-Modification (O15294):** 11 triangulated sites, providing cross-study evidence of extensive auto-regulatory O-GlcNAcylation.
- **Transcriptional Orchestrators:** Classical regulators including **TAB1** (9 sites) and **SP1** (8 sites).

### 3.2 Tissue Specificity in BAP1 Knockout Glycomics (PXD035902)
Evaluating tissue-resolved glycoproteomics in mouse brain and liver reveals marked compartment selectivity:
- **8,514 unique sites** detected in brain versus **2,140 sites** in liver.
- Across 811 sites quantified in both tissues, **61.8% show at least two-fold higher intensity in brain** (median brain/liver ratio: 3.01).
- Site **Thr434 on Q9Z2D6** exhibited **285.7-fold higher intensity** in brain tissue relative to liver.

### 3.3 Quantitative SILAC Dynamics (PXD039536)
Tracking paired Light/Heavy isotope channels in metabolic labeling:
- **393 quantified sites** demonstrated widespread heavy-channel induction (median fold-change 1.32x, 76.1% heavy-biased).
- **Ser1113 on NUP153 (P49790)** was the top dynamic responder (log2 fold-change = 7.39).
- **82 of these responsive SILAC sites** are triangulated across independent external studies.

### 3.4 Cross-Laboratory Quantitative Concordance
Empirical evaluation of cross-study mass spectrometry reproducibility:
- **Intra-laboratory replicates** (e.g. GlycoID serum vs. cytosol replicates from the same publication) correlate strongly ($r = 0.953$, $n = 44$).
- **Biologically related OGT-network studies** show moderate concordance ($r = 0.506$, $n = 88$).
- **Cross-laboratory uncalibrated intensities** between independent facilities (US vs. China cohorts) show weak correlation ($r = 0.274$, $n = 84$).
- **Scientific Takeaway:** Biological site identities harmonize reliably across search engines and laboratories, but uncalibrated mass spectrometry intensity metrics require within-study normalization and cannot be pooled directly without calibration.

---

## 4. Querying and Exploring the Data

### Analytical SQL Queries (DuckDB)
The compendium is stored in columnar Parquet format, enabling instant multi-study queries:

```sql
-- Query all public observations, tissues, and instruments for a specific protein site
SELECT dataset_id, cond_tissue, cond_treatment, metric_name, metric_value, inst_model, prov_country
FROM read_parquet('megatensor/pride/staging/observations.parquet')
WHERE protein_id_raw = 'P51610' AND residue_pos_raw = 579;
```

A complete analytical query suite is provided in [`queries/queries.sql`](queries/queries.sql).

### Machine Learning Feature Matrices
Pre-computed matrices ready for downstream modeling and statistical analysis are stored in [`exports/`](exports/):
- `site_x_condition.parquet`: Normalized site intensity across tissues, treatments, and SILAC channels.
- `site_x_features.parquet`: Site features including sequence windows (+/- 7 residues), protein domain annotations, and evidence scores.

---

## 5. Technical Report and Analytical Findings

The comprehensive undergraduate thesis report: featuring the 7-axis ontology model, cross-laboratory concordance analysis, and 15 embedded analytical figures: is available as a publication-grade PDF:
- **Final Report (PDF):** [`7-Axis-Megatensor-UROP-Final-Report.pdf`](7-Axis-Megatensor-UROP-Final-Report.pdf) (8 pages, fully formatted with embedded vector figures and complete data tables).

---

## 6. Pipeline Reproduction

The complete pipeline can be executed via `just`:

```bash
# Setup environment and dependencies
just setup

# Ingest canonical reference libraries (O-GlcNAc Database & O-GlcNAcAtlas)
just download
just canon

# Ingest curated PRIDE experimental result tables
just pride-discover
just pride-download
just pride-tensorize

# Assemble the 7-axis megatensor and generate analytical results
just union
just analyze
just enrich
just export
```

---

## 7. Citations and Data Acknowledgements

- **The O-GlcNAc Database (MCW):** Wulff-Fuentes et al. *Scientific Data* 8, 25 (2021). [PMID: 33479245](https://pubmed.ncbi.nlm.nih.gov/33479245).
- **O-GlcNAcAtlas 4.0:** Ma et al. *Glycobiology* 31, 719-723 (2021); Hou et al. *Journal of Molecular Biology* 437, 169033 (2025).
- **PRIDE Archive:** European Bioinformatics Institute (EMBL-EBI). Ingested via [`pride-ingest`](https://github.com/filiprumenovski/pride-ingest).
- **Support:** This work was supported by the Wayne State University Undergraduate Research Opportunities Program (UROP) in the Department of Chemistry.
