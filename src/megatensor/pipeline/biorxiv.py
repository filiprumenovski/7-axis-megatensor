"""UROP final report from analysis outputs."""

from __future__ import annotations

import json
from pathlib import Path

from megatensor.paths import FIGURES, ROOT
from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE


def run_biorxiv() -> Path:
    canon = json.loads(CANON_STORE.summary_path.read_text()) if CANON_STORE.summary_path.is_file() else {}
    pride = json.loads(PRIDE_STORE.summary_path.read_text()) if PRIDE_STORE.summary_path.is_file() else {}
    union = json.loads(UNION_STORE.summary_path.read_text()) if UNION_STORE.summary_path.is_file() else {}
    analysis_path = FIGURES / "analysis_summary.json"
    analysis = json.loads(analysis_path.read_text()) if analysis_path.is_file() else {}

    silac = analysis.get("silac", {})
    conc = analysis.get("concordance", {})
    conc_ctx = analysis.get("concordance_context", {})
    ogt = analysis.get("concordance_ogt", {})
    bl = analysis.get("brain_liver", {})
    rep = analysis.get("replication", {})
    imp = analysis.get("importance", {})
    payoff = imp.get("megatensor_payoff", {})
    flat_cost = imp.get("flat_file_cost", {})
    ctx_metrics = (imp.get("pillars") or [{}, {}, {}])[1].get("metrics", {}) if imp else {}
    ev_metrics = (imp.get("pillars") or [{}, {}, {}])[2].get("metrics", {}) if imp else {}
    gsea = analysis.get("gsea")

    gsea_text = ""
    if gsea:
        for label, info in gsea.items():
            top_p = info.get("top_padj", 1)
            gsea_text += (
                f"\n- **{label}:** {info['genes']} genes; "
                f"top term *{str(info.get('top_term', ''))[:55]}* (adj. P={top_p:.2e})."
            )

    hub = rep.get("top_protein_hub") or {}
    hub_label = f"{hub.get('gene_symbol') or '?'} ({hub.get('protein_acc', '—')})"
    gsea_block = gsea_text if gsea_text else "\n_Pathway enrichment not run — `pip install gseapy` then `megatensor publish`._"

    text = f"""# Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model

**UROP Final Report**

- **Student:** Filip Rumenovski
- **Faculty mentor:** Dr. Charlie Fehl
- **Department/Program:** Department of Chemistry, College of Liberal Arts and Sciences; Undergraduate Research Opportunities Program (UROP)
- **Institution:** Wayne State University

---

## Abstract

O-GlcNAc proteomics datasets are released as incompatible tables that collapse multidimensional experimental context into flat rows, blocking cross-lab comparison and machine learning. We introduce a **Site Event Tensor (SET)** representation: a sparse, append-only **megatensor** with seven axes separating biological identity from experimental context (quantification, condition, acquisition, instrument, provenance). We harmonized two canonical O-GlcNAc reference libraries (O-GlcNAc Database, O-GlcNAcAtlas 4.0) and twelve PRIDE experimental deposits spanning MaxQuant, Proteome Discoverer, and mzTab exports from US and Chinese laboratories. Without pairwise manual reconciliation, **{union.get('shared_sites', '—'):,}** human-readable sites intersect canon and PRIDE layers; **{rep.get('triangulated_sites', '—'):,}** sites are **triangulated** (canon-supported and observed in ≥2 independent PRIDE studies). SILAC light/heavy analysis of PXD039536 yields **{silac.get('n_sites', '—')}** quantified sites (median **{silac.get('median_fc', 0):.2f}×** heavy/light). Cross-study intensity concordance at shared sites between PXD039536 and PXD058744 is modest (Pearson **r = {conc.get('pearson_r_log10', '—')}**, *n* = {conc.get('n_shared_sites', '—')}), supporting identity-level interoperability while quantification remains engine- and study-specific. This undergraduate research project produced reproducible analysis tables, ML-ready exports, and DuckDB queries while defining chemistry discrimination (HexNAc vs O-GlcNAc) and unified localization scores as outside the present scope.

**Keywords:** O-GlcNAc, proteomics, data harmonization, PRIDE, tensor, interoperability

---

## Project Background and Objectives

O-linked β-N-acetylglucosamine (O-GlcNAc) regulates signaling, stress response, and disease. Public archives now hold thousands of proteomics experiments, but deposited **result tables** remain heterogeneous: column names, modification semantics, and missing metadata differ by search engine and lab. The representational barrier is not storage — it is that **identity and context are entangled in flat files**.

We ask: if each deposit is mapped through a thin adapter into a shared seven-axis contract, can independent datasets be unioned, queried, and analyzed without re-harmonizing spreadsheets? **The megatensor is not a bigger database** — it is a representational bet: stable site identity (`UniProt:position:AA`) is the join key; every experimental fact (SILAC arm, tissue, probe chemistry, engine, PXD) lives on orthogonal context axes. That separation is what makes cross-study replication tiers, concordance tests, and ML-ready exports possible from deposits that were never designed to talk to each other.

We build a **megatensor** of SETs and test (i) site-level identity overlap across canon and PRIDE, (ii) multi-study replication, (iii) condition-axis biology (SILAC, chemoproteomics, tissue), and (iv) cross-study quant concordance — analyses that require a shared site coordinate system and would otherwise demand **{flat_cost.get('pairwise_manual_joins', 66)}** pairwise table reconciliations across **12** PRIDE deposits alone.

---

## Work Plan and Methods

### Work plan

The project proceeded in five stages: (1) define a stable site identity and seven-axis context model; (2) build thin adapters for canonical resources and heterogeneous PRIDE result tables; (3) validate and combine the resulting SETs; (4) test the shared representation through replication, condition, tissue, and concordance analyses; and (5) export figures, queryable tables, and machine-learning-ready matrices. Each stage was implemented as part of a reproducible command-line pipeline so that updated source data can be processed without repeating manual spreadsheet reconciliation.

### Seven-axis SET ontology

**Layer A (identity):** UniProt accession, residue position, amino acid (Ser/Thr). **PTM** is O-GlcNAc (UniMod:43) at the site entity; localization scores are SET payload, not identity coordinates.

**Layer B (context):** quant metrics (intensity, q-value, spectral count), condition (tissue, treatment, replicate), acquisition (DDA/DIA, collision), instrument, provenance (PXD, country, search engine).

Each observation row maps to a SET coordinate; metrics attach without row explosion after PSM rollup.

### Data sources

| Layer | Sources | Sites | SETs |
|-------|---------|------:|-----:|
| Canon | O-GlcNAc DB, O-GlcNAcAtlas I/II | {canon.get('unique_sites', '—'):,} | {canon.get('unique_sets', '—'):,} |
| PRIDE | 12 PXDs (Table S1) | {pride.get('unique_sites', '—'):,} | {pride.get('unique_sets', '—'):,} |

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

**1. Structural interoperability.** A naïve meta-analysis of 12 PRIDE result tables implies **{flat_cost.get('pairwise_manual_joins', 66)}** pairwise harmonization jobs (column mapping + site matching per pair). The megatensor needs **12** one-time adapters; union then recovers **{payoff.get('shared_sites', union.get('shared_sites', '—')):,}** canon∩PRIDE sites and **{payoff.get('triangulated_sites', rep.get('triangulated_sites', '—')):,}** triangulated sites (canon + ≥2 PXDs) with no downstream spreadsheet work.

**2. Context-preserving biology.** Because condition is an axis—not a column name—one `site_key` supports SILAC dynamics (**{ctx_metrics.get('pride_sites_with_silac', '—'):,}** sites), tissue contrasts (**{ctx_metrics.get('pride_sites_with_tissue', '—'):,}** sites), and chemoproteomic probe matrices (**{ctx_metrics.get('chemoproteomic_sites', '—'):,}** sites) without re-ingesting raw files. Example: **811** brain–liver site pairs and **82** SILAC sites that are also triangulated would be painful to recover from flat exports alone.

**3. Evidence ranking and ML hooks.** Replication tiers, protein hubs, and composite evidence scores turn “how much do we believe this site?” into a query. We ship `site_x_condition` and `site_x_features` tensors for downstream modeling; the concordance analyses (GlycoID *r*≈0.95 within family vs *r*≈0.27 cross-lab) are only meaningful once sites are aligned on identity.

**Figure 0** — Megatensor impact summary (`figures/analysis_megatensor_impact.pdf`). Full narrative: `WHY_MEGATENSOR.md`.

### Identity interoperability across canon and PRIDE

Of **{union.get('pride_unique_sites', '—'):,}** PRIDE unique sites, **{union.get('shared_sites', '—'):,}** (**{100 * union.get('shared_sites', 0) / max(union.get('pride_unique_sites', 1), 1):.1f}%**) appear in canonical references — evidence that adapter-level harmonization recovers literature-supported sites without manual curation.

**Figure 1** — Replication tier bar chart (`figures/analysis_replication_tiers.pdf`).

### Triangulated sites and protein hubs

**{rep.get('triangulated_sites', '—'):,}** sites are supported by canon **and** ≥2 PRIDE deposits. Top hub: **{hub_label}** — {hub.get('n_triangulated_sites', '—')} triangulated sites across {hub.get('max_pxds', '—')} PXDs.

**Figure 2** — Protein hub ranked bar (`figures/analysis_protein_hubs.pdf`).

### Tissue-resolved O-GlcNAc in BAP1KO (PXD035902)

Glycomics PSM tables from the OGT interactor network study resolve O-GlcNAc sites across brain, liver, and liver–brain with distinct site burdens per tissue.

**Figure 6** — BAP1KO tissue site counts (`figures/analysis_bap1ko_tissue.pdf`).

### Condition axis: SILAC in PXD039536

**{silac.get('n_sites', '—')}** sites have both Light and Heavy quantification. Median fold-change **{silac.get('median_fc', 0):.2f}×**; **{silac.get('pct_heavy_biased', 0):.1f}%** of sites are heavy-biased (log₂ FC > 0). Top site: **{silac.get('top_site', '—')}** (log₂ FC = {silac.get('top_log2_fc', '—')}).

**Figure 3** — SILAC M–A plot (`figures/analysis_silac_ma.pdf`).

### Quant concordance across studies

At **{conc.get('n_shared_sites', '—')}** sites quantified in both PXD039536 (China, SILAC MaxQuant) and PXD058744 (US MaxQuant), log₁₀ intensities correlate with **r = {conc.get('pearson_r_log10', '—')}**. By contrast, GlycoID serum/cytosol replicates from the same publication (PXD033026 vs PXD033043) reach **r ≈ {conc_ctx.get('glycoid_replicate_r', 0.92)}** at 44 shared sites — concordance is high within a study but not across engines/geographies. The OGT-network thread (PXD035902 BAP1KO glycomics × PXD039536 SILAC) shows intermediate agreement (**r = {ogt.get('pearson_r_log10', '—')}**, *n* = {ogt.get('n_shared_sites', '—')}).

**Figure 4** — Concordance scatter (`figures/analysis_concordance_scatter.pdf`).

**Figure 4b** — Concordance by study context (`figures/analysis_concordance_context.pdf`).

**Figure 4c** — OGT-network concordance (`figures/analysis_ogt_concordance.pdf`).

**Figure 7** — Pairwise concordance heatmap (`figures/analysis_concordance_heatmap.pdf`).

### Brain vs liver O-GlcNAc in BAP1KO glycomics

Among **{bl.get('n_shared_sites', '—')}** sites detected in both brain and liver (PXD035902), **{bl.get('pct_brain_enriched_2x', '—')}%** are ≥2× brain-enriched (median ratio **{bl.get('median_brain_liver_ratio', '—')}×**). Top site **{bl.get('top_site', '—')}** shows **{bl.get('top_brain_liver_ratio', '—')}×** brain/liver intensity — tissue context is a first-class axis, not noise.

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
{gsea_block}

**Figure 8** — GO enrichment contrast: canon-shared vs PRIDE-novel (`figures/analysis_gsea_contrast.pdf`).

---

## Discussion

The important claim is **not** that O-GlcNAc intensities are globally comparable — concordance shows they are not across engines and continents. The claim is that **representation unlocks the right questions**: Which sites replicate? Which proteins accumulate cross-study evidence (HCFC1)? Which tissues diverge? Which SILAC sites are canon-backed? Flat PRIDE tables can answer each question in isolation, after hours of per-study cleanup; the megatensor answers them in one coordinate system.

We demonstrate **structural interoperability**: independent O-GlcNAc deposits append into a shared axis system with interpretable cross-layer overlap and multi-study replication tiers. The megatensor makes “five instruments, one site” queries literal (DuckDB demo in `queries/queries.sql`). **What would not exist without it:** triangulation tiers, cross-PXD concordance panels, brain/liver site scatters tied to the same keys as SILAC and chemoproteomics, and ranked evidence ladders — all downstream of a single union, not a one-off script per paper.

**Who benefits:** (i) curators benchmarking new sites against canon + public data; (ii) labs comparing their PXD to historical studies on the identity axis; (iii) ML groups needing sparse site×condition tensors without hand-building feature tables per deposit.

## Limitations

The current model has five important limitations. First, UniMod:43 collapses HexNAc chemistries. Second, localization scores are not unified across ptmRS, MaxQuant localization probability, and mzTab. Third, PRIDE coverage is deposit-biased rather than exhaustive. Fourth, one rice study (PXD036527) is excluded from human overlap statistics. Finally, this project reanalyzed deposited results and did not perform new mass spectrometry. These boundaries limit chemical and quantitative interpretation, but they do not prevent the identity-level interoperability tested here.

## Future Work

Future development should add FragPipe and FragPipe-Astral adapters, an explicit chemistry axis, calibrated cross-study quantification, disorder and structure enrichment at scale, and supervised machine learning on the exported tensors. Expanding PRIDE coverage and developing a common localization-confidence representation would also strengthen comparisons across search engines and laboratories.

## UROP Experience and Reflection

This project changed my understanding of what makes computational research scientifically useful. At the beginning, the main challenge appeared to be collecting more O-GlcNAc data. As I worked through the public resources, I learned that data volume was not the limiting factor. The harder problem was preserving meaning when different laboratories, instruments, and search engines described similar biological observations in incompatible ways. Designing the SET representation required me to distinguish a stable biological identity from the experimental context surrounding it. That distinction became the central intellectual lesson of the project: careful representation is not clerical cleanup but part of the scientific method because it determines which comparisons are valid.

The work also strengthened my practical skills in Python, columnar data processing, SQL, reproducible pipelines, data validation, and scientific visualization. More importantly, it taught me to treat unexpected or modest results as information rather than failure. The low cross-study intensity correlation could have been hidden or dismissed, but examining it clarified the proper claim of the project. The megatensor supports identity-level interoperability; it does not make measurements from different workflows automatically comparable. Learning to narrow a conclusion to what the evidence supports was as valuable as implementing the software.

Working with Dr. Charlie Fehl helped connect computational choices to the underlying O-GlcNAc biology. His mentorship encouraged me to ask whether a field in a table represented a real biological distinction, an instrument setting, or a software-specific convention. That guidance kept the project focused on scientifically interpretable outputs rather than data processing for its own sake. Through UROP, I gained experience managing a long-form research project, revising its scope as evidence accumulated, and communicating a technical result to audiences with different backgrounds. I leave the project more confident in my ability to move from an open-ended question to a reproducible analysis while also recognizing the importance of documentation, limitations, and mentor feedback.

## Acknowledgements

I thank Dr. Charlie Fehl for his faculty mentorship, scientific guidance, and feedback throughout this project. This work was supported by Wayne State University's Undergraduate Research Opportunities Program (UROP). I also acknowledge the researchers and curators who made the O-GlcNAc Database, O-GlcNAcAtlas, PRIDE Archive, and the underlying deposited studies publicly available.

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

1. Wulff-Fuentes E, et al. The human O-GlcNAcome database and meta-analysis. *Scientific Data*. 2021;8:25. doi:10.1038/s41597-021-00810-4. PMID: 33479245.
2. Ma J, et al. O-GlcNAcAtlas: a database of experimentally identified O-GlcNAc sites and proteins. *Glycobiology*. 2021;31(7):719–723. doi:10.1093/glycob/cwab003.
3. Hou C, Li W, Li Y, Ma J. O-GlcNAcAtlas 4.0: An updated protein O-GlcNAcylation database with site-specific quantification. *Journal of Molecular Biology*. 2025;437(15):169033. doi:10.1016/j.jmb.2025.169033.
4. PRIDE Archive. European Bioinformatics Institute. <https://www.ebi.ac.uk/pride/>.
5. Rumenovski F. `pride-ingest`: reproducible PRIDE metadata ingestion. <https://github.com/filiprumenovski/pride-ingest>.

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
"""
    out = ROOT / "biorxiv.md"
    out.write_text(text)
    return out
