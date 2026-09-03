# Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model

**UROP Final Report**

- **Student:** Filip Rumenovski
- **Faculty mentor:** Dr. Charlie Fehl
- **Department/Program:** Department of Chemistry, College of Liberal Arts and Sciences; Undergraduate Research Opportunities Program (UROP)
- **Institution:** Wayne State University

---

## Abstract

Public O-GlcNAc proteomics datasets are deposited in incompatible tabular formats. Different search engines and laboratories report conflicting column names, modification notations, and metadata. This fragmentation prevents cross-study comparisons and machine learning. We developed a Site Event Tensor (SET) data model to resolve this issue. The model uses a sparse, append-only megatensor with seven orthogonal axes that separate biological site identity from experimental context. We harmonized two canonical databases (O-GlcNAc Database, O-GlcNAcAtlas 4.0) and twelve PRIDE experimental datasets across MaxQuant, Proteome Discoverer, and mzTab pipelines. Without manual spreadsheet editing, 4,376 unique sites intersect canonical and PRIDE layers. In addition, 353 sites are triangulated, meaning they appear in canonical references and at least two independent PRIDE projects. SILAC quantification of PXD039536 identified 393 quantified sites with a median heavy/light ratio of 1.32. Intensity concordance between independent projects (PXD039536 and PXD058744) is low (Pearson r = 0.274, n = 84). This demonstrates that biological site identity harmonizes across studies, while raw intensity scales remain specific to individual laboratories and instruments. The pipeline outputs reproducible Parquet tables, DuckDB queries, and machine learning tensors.

**Keywords:** O-GlcNAc, proteomics, data harmonization, PRIDE, tensor, interoperability

---

## 1. Project Background and Objectives

O-linked beta-N-acetylglucosamine (O-GlcNAc) is a dynamic post-translational modification that controls cell signaling, transcription, and metabolic stress responses. Public repositories store thousands of mass spectrometry runs, but published result tables remain fragmented. Different search engines (such as MaxQuant, Proteome Discoverer, and DIA-NN) output different table schemas. PTM site localization probabilities and peptide definitions are stored inconsistently. In standard flat files, biological site identity and experimental conditions are merged into single text strings.

To compare datasets, researchers typically write custom scripts for each pair of publications. For 12 datasets, pairwise comparison requires 66 separate conversion steps.

This project tests whether a unified seven-axis data contract can eliminate manual reconciliation. The core concept is simple: biological site identity (`UniProt:position:AA`) serves as a single join key. All experimental variables (tissue, chemical probe, SILAC label, instrument model, search engine, and PRIDE accession) are assigned to independent context axes. This structure enables three direct tests:
1. Measuring site overlap between canonical reference databases and raw PRIDE deposits.
2. Identifying high-confidence replicated sites across independent laboratories.
3. Testing whether quantitative intensity measurements correlate across different mass spectrometry facilities.

---

## 2. Work Plan and Methods

### 2.1 Work Plan

The project was executed in five technical stages:
1. Define the seven-axis schema separating site identity from experimental variables.
2. Build automated ingestion adapters for canonical databases and PRIDE result tables.
3. Validate and combine individual datasets into a unified Parquet store.
4. Run cross-study analyses covering replication tiers, SILAC ratios, tissue contrasts, and quantitative concordance.
5. Export structured analytical tables, database views, and figures.

All steps were integrated into an automated command-line workflow using a Makefile/Justfile.

### 2.2 Seven-Axis SET Schema

The data model divides each observation into two layers:
- **Layer A (Site Identity):** UniProt protein accession, residue position, and amino acid (Ser or Thr). The chemical modification is fixed as O-GlcNAc (UniMod:43). PTM localization scores are stored as observation payloads rather than identity keys.
- **Layer B (Experimental Context):** Seven orthogonal context axes: (1) quantification (intensity, signal-to-noise ratio, fold change, q-value); (2) biological condition (tissue type, disease state, treatment, biological replicate); (3) chemical probe (enrichment tag, cleavable linker chemistry); (4) acquisition (instrument mode DDA vs. DIA, collision energy, fragmentation type); (5) instrument model; (6) search engine (identification software, scoring metric, database parameters); and (7) provenance (PRIDE project accession, sample accession, country, publication DOI).

### 2.3 Data Sources

| Layer | Sources | Unique Sites | SET Observations |
|-------|---------|-------------:|-----------------:|
| Canon | O-GlcNAc DB, O-GlcNAcAtlas I/II | 43,853 | 81,486 |
| PRIDE | 12 PXDs (Orbitrap Elite, Fusion, Lumos; MaxQuant, PD, mzTab) | 13,567 | 36,737 |

We integrated two canonical reference databases and twelve PRIDE deposits. We processed deposited identification and quantification tables directly. We did not re-search raw mass spectrometer files. PRIDE projects were chosen to maximize diversity across instrument vendors, search software, and geographic locations.

### 2.4 Software Implementation

The pipeline was implemented in Python 3.11 using Polars for columnar data manipulation and DuckDB for SQL analytics. The pipeline executes sequentially:
`megatensor canon` -> `pride-tensorize` -> `union` -> `analyze` -> `enrich` -> `export`.
Source code and documentation are hosted at: <https://github.com/filiprumenovski/7-axis-megatensor>.

---

## 3. Results

### 3.1 Engineering Payoffs of the Megatensor Architecture

The seven-axis structure provides three distinct operational advantages:
1. **Reduced integration effort:** Comparing 12 PRIDE deposits pairwise requires 66 individual conversion scripts. In the megatensor framework, each repository requires only one ingestion adapter (12 adapters total). A single union operation identified 4,376 shared canonical-PRIDE sites and 353 triangulated sites without manual table formatting.
2. **Context preservation:** Biological conditions remain queryable dimensions rather than unstructured column labels. A single site key connects 1,035 SILAC sites, 8,514 tissue-resolved sites, and 4,622 chemoproteomic sites. For instance, the pipeline recovered 811 brain-liver site pairs and 82 triangulated SILAC sites through basic SQL queries.
3. **Structured evidence ranking:** Sites can be ranked directly by replication depth, hub connectivity, and spectral evidence. The system exports clean feature matrices (`site_x_condition.parquet` and `site_x_features.parquet`) for machine learning.

**Figure 0**: Megatensor impact summary (`figures/analysis_megatensor_impact.pdf`).

### 3.2 Cross-Layer Site Recovery and Replication Tiers

Of the 13,567 unique sites extracted from PRIDE deposits, 4,376 (32.3%) match canonical literature references. This confirms that automated ingestion adapters extract valid biological sites without manual curation.
We classified all identified sites into four confidence tiers:
- **PRIDE-only:** 9,191 sites observed in only one PRIDE project.
- **PRIDE multi-study:** 1,211 sites confirmed by two or more PRIDE accessions.
- **Canon-intersected:** 4,376 sites present in both PRIDE and canonical databases.
- **Triangulated:** 353 high-confidence sites supported by canonical literature and at least two independent PRIDE deposits.

**Figure 1**: Replication tiers (`figures/analysis_replication_tiers.pdf`).

### 3.3 High-Confidence Protein Hubs

The 353 triangulated sites cluster on specific regulatory proteins. HCFC1 (Host Cell Factor 1, P51610) is the top hub, containing 18 triangulated sites confirmed across 3 independent PRIDE studies. Other prominent hubs include nuclear pore complex proteins NUP214 (14 sites) and NUP98 (12 sites), O-GlcNAc transferase (OGT, 11 sites), TAB1 (9 sites), and transcription factor SP1 (8 sites).

**Figure 2**: Top protein hubs (`figures/analysis_protein_hubs.pdf`).

### 3.4 Tissue-Specific O-GlcNAcylation in BAP1 Knockout (PXD035902)

Analysis of PXD035902 resolved 8,514 unique sites in mouse brain and 2,140 sites in liver.
A subset of 811 sites was quantified in both tissues. Of these shared sites, 61.8% exhibited at least two-fold higher intensity in brain (median brain/liver ratio: 3.01). The largest difference was observed on Q9Z2D6 at Thr434, which showed a 285.7-fold higher intensity in brain tissue. This confirms that tissue differences reflect real biological regulation rather than baseline technical variation.

**Figure 6**: BAP1KO tissue burden (`figures/analysis_bap1ko_tissue.pdf`).  
**Figure 6b**: Brain vs. liver contrast (`figures/analysis_bap1ko_brain_liver.pdf`).

### 3.5 SILAC Dynamics (PXD039536)

In PXD039536, 393 sites had paired light and heavy isotope measurements. The median fold change was 1.32, with 76.1% of quantified sites exhibiting positive log2 fold changes (heavy-biased). The top responder was P49790 at Ser1113 (log2 FC = 7.39). Cross-referencing against the full megatensor revealed that 82 of these 393 SILAC sites are also triangulated across canonical databases and additional PRIDE studies.

**Figure 3**: SILAC M-A plot (`figures/analysis_silac_ma.pdf`).  
**Figure 3b**: SILAC triangulation overlay (`figures/analysis_silac_triangulation.pdf`).

### 3.6 Cross-Study Quantitative Concordance

We tested whether raw mass spectrometry intensities can be compared directly across independent laboratories.
At 84 shared sites quantified in both PXD039536 (China, SILAC MaxQuant) and PXD058744 (United States, label-free MaxQuant), log10 intensities showed weak correlation (r = 0.274).
In contrast, replicates from the same publication and laboratory (PXD033026 vs. PXD033043, GlycoID serum and cytosol) showed strong concordance (r = 0.953, n = 44).
Studies examining the OGT interaction network (PXD035902 and PXD039536) showed moderate concordance (r = 0.506, n = 88).
These results demonstrate that site identities align reliably across studies, but uncalibrated intensity values cannot be compared across different laboratories.

**Figure 4**: Concordance scatter (`figures/analysis_concordance_scatter.pdf`).  
**Figure 4b**: Concordance by context (`figures/analysis_concordance_context.pdf`).  
**Figure 4c**: OGT-network concordance (`figures/analysis_ogt_concordance.pdf`).  
**Figure 7**: Pairwise concordance heatmap (`figures/analysis_concordance_heatmap.pdf`).

### 3.7 Evidence Scores and Chemoproteomic Probe Partitioning

We computed composite evidence scores based on canonical database presence, PRIDE study replication count, and spectral observation frequency. The highest-scoring sites were P51610:579:T (HCFC1), Q14119:621:S (BAP1 cofactor), and P49790:1113:S (NUP153).
In PXD063995, three distinct chemical probes (PC, DDE, and DADPS) were used to enrich modified peptides. Individual O-GlcNAc sites showed strong probe preferences, demonstrating that chemical enrichment methods introduce distinct labeling biases.

**Figure 9**: Evidence ladder (`figures/analysis_evidence_ladder.pdf`).  
**Figure 5**: Chemoproteomics specificity (`figures/analysis_chemoproteomics.pdf`).

### 3.8 Pathway Enrichment and Triangulated Intensity Panel

Figure 10 shows intensity measurements for top triangulated sites across individual PRIDE projects. Gene Ontology enrichment (Enrichr GO Biological Process) identified positive regulation of DNA-templated transcription as the top enriched term in both canonical sites (p_adj = 1.13e-16) and PRIDE-novel sites (p_adj = 8.15e-08). Novel sites were also enriched for cellular stress responses, protein folding, and intracellular transport.

**Figure 10**: Triangulated site panel (`figures/analysis_triangulated_heatmap.pdf`).  
**Figure 8**: GO pathway enrichment (`figures/analysis_gsea_contrast.pdf`).

---

## 4. Discussion

This project demonstrates that biological site identity can be harmonized across disparate proteomics repositories without manual table editing. The megatensor framework converts heterogeneous flat files into a queryable relational store.

The concordance analysis provides an important practical conclusion: while site identities are consistent across studies, raw quantitative values are not directly comparable across different laboratories and instruments. Variations in sample preparation, instrument tuning, and search engine algorithms produce substantial baseline shifts. Consequently, meta-analyses should focus on identity overlap, replication counts, and within-experiment relative fold changes rather than raw pooled intensities.

The resulting database directly supports three use cases:
1. **Benchmarking:** Curators can quickly test whether newly observed sites have prior literature or public repository support.
2. **Cross-study verification:** Experimentalists can determine if an identified site replicates across external datasets.
3. **Machine learning:** Computational groups can train predictive models on standardized site-by-feature matrices.

---

## 5. Limitations

The current implementation has five technical boundaries:
1. UniMod:43 does not distinguish O-GlcNAc from isobaric HexNAc stereoisomers such as O-GalNAc.
2. Localization confidence scores (ptmRS, MaxQuant localization probability, and mzTab scores) are not yet converted into a single unified probability scale.
3. PRIDE coverage was restricted to twelve targeted deposits rather than the full public archive.
4. One rice dataset (PXD036527) was excluded from human site overlap statistics.
5. All analyses were performed on published result tables without re-searching raw mass spectrometry data files.

---

## 6. Future Work

Planned technical improvements include:
1. Developing automated ingestion adapters for FragPipe and FragPipe-Astral pipelines.
2. Adding a dedicated chemistry axis to track isobaric glycan definitions.
3. Implementing cross-study intensity normalization algorithms.
4. Integrating protein structure and disorder annotations from AlphaFold.
5. Training machine learning classifiers on the exported feature matrices.

---

## 7. UROP Experience and Reflection

This UROP project provided practical experience in scientific computing, data engineering, and mass spectrometry proteomics.

At the beginning of the internship, I expected data collection to be the primary challenge. Working with public archives quickly demonstrated that data harmonization was the real difficulty. Public proteomics results are published in conflicting table layouts with different column names and modification notations. Designing the SET schema taught me how to separate permanent biological identifiers from variable experimental conditions.

A key practical lesson was interpreting low correlation values as meaningful scientific findings rather than software defects. The low intensity correlation across laboratories (r = 0.274) clarified the actual capabilities of the pipeline: the system provides reliable identity harmonization, but raw mass spectrometry intensities require study-specific calibration.

Under the mentorship of Dr. Charlie Fehl, I learned how software engineering choices connect to chemical and biological reality. Dr. Fehl helped me distinguish meaningful biological differences from instrument settings and software artifacts. Through this project, I gained substantial experience in Python, Polars, DuckDB, Parquet storage, automated pipelines, and technical documentation.

---

## Acknowledgements

I thank Dr. Charlie Fehl for his mentorship, technical guidance, and feedback throughout this project. This work was supported by Wayne State University's Undergraduate Research Opportunities Program (UROP). I also thank the authors and maintainers of the O-GlcNAc Database, O-GlcNAcAtlas, and PRIDE Archive for making their data publicly accessible.

---

## Data Availability

- **Canonical Downloads:** O-GlcNAc Database, O-GlcNAcAtlas 4.0.
- **PRIDE Accessions:** Complete accession list in `figures/pride_glyco_picks.csv`.
- **Analysis Tables:** Processed Parquet tables in `megatensor/analysis/*.parquet`.
- **ML Tensors:** Exported to `exports/site_x_condition.parquet` and `exports/site_x_features.parquet`.

## Code Availability

Pipeline and adapters in this repository. Reproduce via:
```bash
just setup && just download && just canon
just pride-download && just pride-tensorize
just union && just analyze && just enrich && just export
```

---

## References

1. Wulff-Fuentes E, et al. The human O-GlcNAcome database and meta-analysis. *Scientific Data*. 2021;8:25. doi:10.1038/s41597-021-00810-4.
2. Ma J, et al. O-GlcNAcAtlas: a database of experimentally identified O-GlcNAc sites and proteins. *Glycobiology*. 2021;31(7):719–723. doi:10.1093/glycob/cwab003.
3. Hou C, Li W, Li Y, Ma J. O-GlcNAcAtlas 4.0: An updated protein O-GlcNAcylation database with site-specific quantification. *Journal of Molecular Biology*. 2025;437(15):169033. doi:10.1016/j.jmb.2025.169033.
4. PRIDE Archive. European Bioinformatics Institute. <https://www.ebi.ac.uk/pride/>.
5. Rumenovski F. *pride-ingest*: Reproducible PRIDE metadata ingestion. Available at: <https://github.com/filiprumenovski/pride-ingest>.

---

## Figure List

| Fig | File | Key Scientific Finding |
|-----|------|------------------------|
| 0 | `analysis_megatensor_impact` | Architecture reduction (66 manual joins -> 12 adapters) and multi-axis preservation. |
| 1 | `analysis_replication_tiers` | Stratification of 13,567 PRIDE sites into replication and canonical intersection tiers. |
| 2 | `analysis_protein_hubs` | Cross-study evidence accumulation on key hubs (HCFC1, NUP214, NUP98, OGT). |
| 3 | `analysis_silac_ma` | PXD039536 SILAC condition axis dynamics (393 sites; 76.1% heavy-biased; median 1.32x). |
| 3b | `analysis_silac_triangulation` | Cross-layer overlay isolating 82 triangulated sites within quantitative SILAC space. |
| 4 | `analysis_concordance_scatter` | Weak quantitative correlation (r = 0.274) between US and China MaxQuant runs at 84 shared sites. |
| 4b | `analysis_concordance_context` | Stratified concordance: intra-lab replicates (r ~ 0.95) vs. cross-lab cohorts (r ~ 0.27). |
| 4c | `analysis_ogt_concordance` | Moderate concordance (r = 0.506) across biologically related OGT-network experiments. |
| 5 | `analysis_chemoproteomics` | Chemoproteomic probe partitioning demonstrating probe-specific site selectivity. |
| 6 | `analysis_bap1ko_tissue` | Total O-GlcNAc site burden partitioned across brain, liver, and shared tissue subsets. |
| 6b | `analysis_bap1ko_brain_liver` | Site-specific brain enrichment (61.8% >= 2x; Q9Z2D6:434:T at 285.7x) across 811 shared sites. |
| 7 | `analysis_concordance_heatmap` | Global cross-study pairwise Pearson correlation matrix across all 12 PRIDE deposits. |
| 8 | `analysis_gsea_contrast` | Functional GO Biological Process contrast between canon-shared and PRIDE-novel gene sets. |
| 9 | `analysis_evidence_ladder` | Multi-layer composite evidence ranking isolating top O-GlcNAc sites for downstream ML. |
| 10 | `analysis_triangulated_heatmap` | Multi-PXD intensity observation matrix across top high-confidence triangulated sites. |
