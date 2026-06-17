# Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model

**Filip Rumenovski¹, Charlie Fehl¹**  
¹ Washington State University

*Draft preprint — generated from reproducible pipeline outputs.*

---

## Abstract

O-GlcNAc proteomics datasets are released as incompatible tables that collapse multidimensional experimental context into flat rows, blocking cross-lab comparison and machine learning. We introduce a **Site Event Tensor (SET)** representation: a sparse, append-only **megatensor** with seven axes separating biological identity from experimental context (quantification, condition, acquisition, instrument, provenance). We harmonized two canonical O-GlcNAc reference libraries (O-GlcNAc Database, O-GlcNAcAtlas 4.0) and twelve PRIDE experimental deposits spanning MaxQuant, Proteome Discoverer, and mzTab exports from US and Chinese laboratories. Without pairwise manual reconciliation, **4,376** human-readable sites intersect canon and PRIDE layers; **353** sites are **triangulated** (canon-supported and observed in ≥2 independent PRIDE studies). SILAC light/heavy analysis of PXD039536 yields **393** quantified sites (median **1.32×** heavy/light). Cross-study intensity concordance at shared sites between PXD039536 and PXD058744 is modest (Pearson **r = 0.274**, *n* = 84), supporting identity-level interoperability while quantification remains engine- and study-specific. We release analysis tables, ML-ready exports, and DuckDB queries. This is a methods-and-resource paper; chemistry discrimination (HexNAc vs O-GlcNAc) and unified localization scores are explicitly out of scope.

**Keywords:** O-GlcNAc, proteomics, data harmonization, PRIDE, tensor, interoperability

---

## Introduction

O-linked β-N-acetylglucosamine (O-GlcNAc) regulates signaling, stress response, and disease. Public archives now hold thousands of proteomics experiments, but deposited **result tables** remain heterogeneous: column names, modification semantics, and missing metadata differ by search engine and lab. The representational barrier is not storage — it is that **identity and context are entangled in flat files**.

We ask: if each deposit is mapped through a thin adapter into a shared seven-axis contract, can independent datasets be unioned, queried, and analyzed without re-harmonizing spreadsheets? **The megatensor is not a bigger database** — it is a representational bet: stable site identity (`UniProt:position:AA`) is the join key; every experimental fact (SILAC arm, tissue, probe chemistry, engine, PXD) lives on orthogonal context axes. That separation is what makes cross-study replication tiers, concordance tests, and ML-ready exports possible from deposits that were never designed to talk to each other.

We build a **megatensor** of SETs and test (i) site-level identity overlap across canon and PRIDE, (ii) multi-study replication, (iii) condition-axis biology (SILAC, chemoproteomics, tissue), and (iv) cross-study quant concordance — analyses that require a shared site coordinate system and would otherwise demand **66** pairwise table reconciliations across **12** PRIDE deposits alone.

---

## Methods

### Seven-axis SET ontology

**Layer A (identity):** UniProt accession, residue position, amino acid (Ser/Thr). **PTM** is O-GlcNAc (UniMod:43) at the site entity; localization scores are SET payload, not identity coordinates.

**Layer B (context):** quant metrics (intensity, q-value, spectral count), condition (tissue, treatment, replicate), acquisition (DDA/DIA, collision), instrument, provenance (PXD, country, search engine).

Each observation row maps to a SET coordinate; metrics attach without row explosion after PSM rollup.

### Data sources

| Layer | Sources | Sites | SETs |
|-------|---------|------:|-----:|
| Canon | O-GlcNAc DB, O-GlcNAcAtlas I/II | 43,853 | 81,486 |
| PRIDE | 12 PXDs (Table S1) | 13,567 | 36,737 |

PRIDE picks prioritized engine, instrument, and geography heterogeneity; only deposited result tables were parsed (no raw file re-search).

### Software

Python 3.11, Polars, DuckDB. Pipeline: `megatensor canon` → `pride-tensorize` → `union` → `analyze` → `enrich` → `export`. Code: [repository].

### Analyses

1. **Replication tiers** — PRIDE-only, multi-PXD, canon∩PRIDE, triangulated (canon + ≥2 PXDs).
2. **SILAC** — PXD039536: per-site Heavy/Light mean intensity → log₂ fold-change; M–A plot.
3. **Concordance** — log₁₀ mean intensity at sites observed in both PXD039536 and PXD058744.
4. **Chemoproteomics** — PXD063995 probe matrix (PC, DDE, DADPS, …).
5. **Pathway enrichment** — Enrichr on gene symbols for PRIDE-novel vs canon-shared sites (optional).

---

## Results

### Why the megatensor matters (three payoffs)

**1. Structural interoperability.** A naïve meta-analysis of 12 PRIDE result tables implies **66** pairwise harmonization jobs (column mapping + site matching per pair). The megatensor needs **12** one-time adapters; union then recovers **4,376** canon∩PRIDE sites and **353** triangulated sites (canon + ≥2 PXDs) with no downstream spreadsheet work.

**2. Context-preserving biology.** Because condition is an axis—not a column name—one `site_key` supports SILAC dynamics (**1,035** sites), tissue contrasts (**8,514** sites), and chemoproteomic probe matrices (**4,622** sites) without re-ingesting raw files. Example: **811** brain–liver site pairs and **82** SILAC sites that are also triangulated would be painful to recover from flat exports alone.

**3. Evidence ranking and ML hooks.** Replication tiers, protein hubs, and composite evidence scores turn “how much do we believe this site?” into a query. We ship `site_x_condition` and `site_x_features` tensors for downstream modeling; the concordance analyses (GlycoID *r*≈0.95 within family vs *r*≈0.27 cross-lab) are only meaningful once sites are aligned on identity.

**Figure 0** — Megatensor impact summary (`figures/analysis_megatensor_impact.pdf`). Full narrative: `WHY_MEGATENSOR.md`.

### Identity interoperability across canon and PRIDE

Of **13,567** PRIDE unique sites, **4,376** (**32.3%**) appear in canonical references — evidence that adapter-level harmonization recovers literature-supported sites without manual curation.

**Figure 1** — Replication tier bar chart (`figures/analysis_replication_tiers.pdf`).

### Triangulated sites and protein hubs

**353** sites are supported by canon **and** ≥2 PRIDE deposits. Top hub: **HCFC1 (P51610)** — 18 triangulated sites across 3 PXDs.

**Figure 2** — Protein hub ranked bar (`figures/analysis_protein_hubs.pdf`).

### Tissue-resolved O-GlcNAc in BAP1KO (PXD035902)

Glycomics PSM tables from the OGT interactor network study resolve O-GlcNAc sites across brain, liver, and liver–brain with distinct site burdens per tissue.

**Figure 6** — BAP1KO tissue site counts (`figures/analysis_bap1ko_tissue.pdf`).

### Condition axis: SILAC in PXD039536

**393** sites have both Light and Heavy quantification. Median fold-change **1.32×**; **76.1%** of sites are heavy-biased (log₂ FC > 0). Top site: **P49790:1113:S** (log₂ FC = 7.39).

**Figure 3** — SILAC M–A plot (`figures/analysis_silac_ma.pdf`).

### Quant concordance across studies

At **84** sites quantified in both PXD039536 (China, SILAC MaxQuant) and PXD058744 (US MaxQuant), log₁₀ intensities correlate with **r = 0.274**. By contrast, GlycoID serum/cytosol replicates from the same publication (PXD033026 vs PXD033043) reach **r ≈ 0.9526666666666666** at 44 shared sites — concordance is high within a study but not across engines/geographies. The OGT-network thread (PXD035902 BAP1KO glycomics × PXD039536 SILAC) shows intermediate agreement (**r = 0.506**, *n* = 88).

**Figure 4** — Concordance scatter (`figures/analysis_concordance_scatter.pdf`).

**Figure 4b** — Concordance by study context (`figures/analysis_concordance_context.pdf`).

**Figure 4c** — OGT-network concordance (`figures/analysis_ogt_concordance.pdf`).

**Figure 7** — Pairwise concordance heatmap (`figures/analysis_concordance_heatmap.pdf`).

### Brain vs liver O-GlcNAc in BAP1KO glycomics

Among **811** sites detected in both brain and liver (PXD035902), **61.8%** are ≥2× brain-enriched (median ratio **3.01×**). Top site **Q9Z2D6:434:T** shows **285.7×** brain/liver intensity — tissue context is a first-class axis, not noise.

**Figure 6b** — Brain vs liver scatter (`figures/analysis_bap1ko_brain_liver.pdf`).

### Highest-evidence sites

Composite scoring (canon depth + PRIDE replication + SET support) ranks **P51610:579:T** (HCFC1) and **Q14119** cofactor sites at the top; **82** SILAC-quantified sites are triangulated (canon + ≥2 PXDs).

**Figure 9** — Evidence ladder (`figures/analysis_evidence_ladder.pdf`).

**Figure 10** — Triangulated site × PXD intensity panel (`figures/analysis_triangulated_heatmap.pdf`).

**Figure 3b** — SILAC with triangulation overlay (`figures/analysis_silac_triangulation.pdf`).

### Chemoproteomic probe partitioning

PXD063995 provides multi-probe O-GlcNAc chemoproteomics; site × probe heatmaps show probe-specific intensity patterns consistent with differential labeling chemistry.

**Figure 5** — Chemoproteomics heatmap (`figures/analysis_chemoproteomics.pdf`).

### Pathway context

- **pride_novel:** 800 genes; top term *Positive Regulation Of DNA-templated Transcription (GO:* (adj. P=8.15e-08).
- **canon_shared:** 800 genes; top term *Positive Regulation Of DNA-templated Transcription (GO:* (adj. P=1.13e-16).

**Figure 8** — GO enrichment contrast: canon-shared vs PRIDE-novel (`figures/analysis_gsea_contrast.pdf`).

---

## Discussion

The important claim is **not** that O-GlcNAc intensities are globally comparable — concordance shows they are not across engines and continents. The claim is that **representation unlocks the right questions**: Which sites replicate? Which proteins accumulate cross-study evidence (HCFC1)? Which tissues diverge? Which SILAC sites are canon-backed? Flat PRIDE tables can answer each question in isolation, after hours of per-study cleanup; the megatensor answers them in one coordinate system.

We demonstrate **structural interoperability**: independent O-GlcNAc deposits append into a shared axis system with interpretable cross-layer overlap and multi-study replication tiers. The megatensor makes “five instruments, one site” queries literal (DuckDB demo in `queries/queries.sql`). **What would not exist without it:** triangulation tiers, cross-PXD concordance panels, brain/liver site scatters tied to the same keys as SILAC and chemoproteomics, and ranked evidence ladders — all downstream of a single union, not a one-off script per paper.

**Who benefits:** (i) curators benchmarking new sites against canon + public data; (ii) labs comparing their PXD to historical studies on the identity axis; (iii) ML groups needing sparse site×condition tensors without hand-building feature tables per deposit.

**Limitations (explicit):** (1) UniMod:43 collapses HexNAc chemistries; (2) localization scores are not unified across ptmRS, MaxQuant loc prob, and mzTab; (3) PRIDE coverage is deposit-biased, not exhaustive; (4) one rice study (PXD036527) is excluded from human overlap stats; (5) no new mass spectrometry was performed.

**Future work:** FragPipe/FragPipe-Astral adapters, chemistry axis, calibrated quant harmonization, disorder/structure enrichment at scale, supervised ML on exported tensors.

---

## Data availability

- Canon downloads: O-GlcNAc Database, O-GlcNAcAtlas  
- PRIDE: PXD accessions in `figures/pride_glyco_picks.csv`  
- Analysis tables: `megatensor/analysis/*.parquet`  
- ML exports: `exports/site_x_condition.parquet`, `exports/site_x_features.parquet`

## Code availability

Pipeline and adapters in this repository. Reproduce:

```bash
just setup && just download && just canon
just pride-download && just pride-tensorize
just union && just analyze && just enrich && just export
.venv/bin/python -c "from megatensor.viz.analysis_plots import render_analysis_figures; render_analysis_figures()"
```

## References

1. Wulff-Fuentes et al. *Sci Data* 8, 25 (2021) — O-GlcNAc Database  
2. Ma et al. *Glycobiology* 31, 719–723 (2021); Hou et al. *J Mol Biol* (2025) — O-GlcNAcAtlas  
3. PRIDE Archive / pride-ingest  

---

## Figure list

| Fig | File | Claim |
|-----|------|-------|
| 0 | `analysis_megatensor_impact` | Why structure matters: cost vs payoff |
| 1 | `analysis_replication_tiers` | Replication tiers quantify interoperability depth |
| 2 | `analysis_protein_hubs` | Some proteins accumulate cross-study O-GlcNAc evidence |
| 3 | `analysis_silac_ma` | SILAC condition axis resolves site-level dynamics |
| 4 | `analysis_concordance_scatter` | Shared IDs ≠ shared quant scale |
| 5 | `analysis_chemoproteomics` | Multi-probe chemistry visible on condition axis |
| 6 | `analysis_bap1ko_tissue` | Tissue axis in BAP1KO network study |
| 7 | `analysis_concordance_heatmap` | Study-pair quant concordance overview |
| 8 | `analysis_gsea_contrast` | Functional contrast novel vs shared sites |
| 9 | `analysis_evidence_ladder` | Composite cross-layer evidence ranking |
| 10 | `analysis_triangulated_heatmap` | Triangulated sites across PXDs |
| 3b | `analysis_silac_triangulation` | SILAC FC with triangulation overlay |
| 4b | `analysis_concordance_context` | Within-pipeline vs cross-lab concordance |
| 4c | `analysis_ogt_concordance` | OGT-network study agreement |
| 6b | `analysis_bap1ko_brain_liver` | Tissue-specific O-GlcNAc in BAP1KO |
