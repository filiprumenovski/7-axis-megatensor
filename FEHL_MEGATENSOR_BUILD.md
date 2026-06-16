# FEHL Megatensor: Autonomous Build Doctrine

**Deliverable:** UROP project, "Toward Interoperable O-GlcNAc Proteomics: A Tiered 7-Axis Megatensor Model" (Filip Rumenovski / Dr. Charlie Fehl, WSU).
**Mode:** autonomous build with sparse human checkpoints. The agent runs the loop unsupervised between checkpoints and only stops at the gates in §12.
**Scope guard:** machine learning is **out of scope for this run**. Produce ML-*ready* exports only (tensor shapes, feature matrices). Do not train, fit, or evaluate any model. The hook stays open for a later phase.
**Style:** dense, directive, reproducible. Everything is Parquet, append-only, sparse. Never densify the full tensor.

---

## 0. Definition of done

At the end of the run, the following exist and are reproducible from one command:

1. A populated Megatensor unioning the two canonical O-GlcNAc reference databases plus 4-5 fresh PRIDE result datasets chosen for maximum heterogeneity (instrument, country, search engine), with no manual harmonization downstream of the per-file adapters.
2. **Figure A (axis completeness):** per-source coverage across all seven axes, plus SET count, protein count, site count, metric count.
3. **Figure B (cross-dataset interoperability):** UpSet of site membership across all sources, with the canon-vs-canon reconciliation called out explicitly.
4. **Figure C (single-site multi-condition trajectory):** one biologically meaningful O-GlcNAc site sliced across conditions, pulled live from the union. Must originate from a conditioned PRIDE dataset (canon cannot supply a trajectory).
5. **Enriched Megatensor:** disorder, domain/region, sequence window, and pathway enrichment attached as SET/identity payloads. Structure-derived SASA scoped to visualized sites only.
6. **Live query demo:** a set of DuckDB queries that slice the union by any axis in one statement (the "five labs, five instruments, one site" query is the headline).
7. **Tensor exports:** dense site x condition matrix and an enriched site x feature matrix written to disk (Parquet + .npy). No modeling.
8. **Technical report (`report.md`)** and a clean repo with `README.md` and a one-command runner.

A canon-only version of items 1-3 and 8 is **already a complete, defensible deliverable**. PRIDE and enrichment are additive. Sequence the build (§12) so a short night still lands something whole.

---

## 1. Project canon (carry-through from the two source documents)

**Problem.** O-GlcNAc proteomics data are released as heterogeneous CSV/Excel tables that collapse multidimensional experimental context into two-dimensional rows. The barrier is representational, not tooling: variation in quantification method, acquisition mode, instrument, and provenance makes datasets incomparable without manual effort, blocking reproducibility, cross-lab validation, and ML.

**Thesis.** Separate stable biological identity from variable experimental context, encode every measurement as a Site Event Tensor (SET) in a tiered 7-axis ontology, and union SETs into a sparse, append-only Megatensor. Interoperability then happens structurally, at the representation level, with no pairwise schema reconciliation.

**Objectives (from proposal).**
1. Build a Megatensor by tensorizing the canonical databases plus 2-3+ external datasets.
2. Demonstrate interoperability via cross-dataset comparisons, UpSet plots, and axis-completeness metrics.
3. Enable biological exploration by extracting multi-condition trajectories for selected sites.
4. Prepare ML-ready output formats (dense vectors, embedding matrices, graph inputs). *Export only this run.*
5. Ship documentation and open-source code suitable for lab adoption.

**Success criteria.** SETs from independent datasets append into one Megatensor without manual harmonization; researchers can extract meaningful biological slices (per protein / site / condition); computed tensors are directly exportable for downstream modeling.

**Honest framing for the defense.** The claim is not "zero manual work." Each source needs a thin column-binding adapter (a few lines mapping that file's columns to the canonical contract). The harmonization *into the 7-axis space* is automatic, and that automatic part is the hard problem nobody else solved. State it this way.

---

## 2. The tiered 7-axis ontology

Two semantic layers. Layer A is universal biological identity (the registry). Layer B is experimental context (the append axes).

**Layer A: global biological identity**

1. **Identity.** Protein, isoform, residue position, amino acid. Anchors all downstream representation and enables inter-dataset joinability. All positions are mapped to UniProt canonical sequence coordinates.
2. **PTM.** Modification class (here O-GlcNAc, i.e. HexNAc on Ser/Thr). Localization confidence is recorded but kept **out of the identity coordinate** (see §9 refinement).

**Layer B: experimental context (append axes)**

3. **Quant.** Numeric measurements (intensity, ratio, score, abundance, quantification event), raw or normalized, multiple metrics per event.
4. **Condition.** Biological/experimental perturbation: cell line, tissue, timepoint, treatment, fraction (sup/pel), replicate type.
5. **Acquisition.** MS method choices: DDA/DIA, MS2/MS3, gradient, collision type, enrichment chemistry (lectin/antibody/chemoenzymatic/metabolic).
6. **Instrument.** Hardware identity: vendor, model, source type, PSI-MS CV accession.
7. **Provenance.** Dataset/run identifiers, lab origin, country, DOI/accession, search software + version, preprocessing hash, tensorization version.

`Identity x PTM` defines the canonical **site entity**. The full coordinate `Identity x PTM x Quant x Condition x Acquisition x Instrument x Provenance` defines a single measurement event.

---

## 3. Site Event Tensor and Megatensor

**SET.** The atomic unit. One SET = one site entity (Identity x PTM) measured under one experimental configuration (Condition x Acquisition x Instrument x Provenance), carrying one or more Quant payloads. Materialized as an `indices` object (the axis IDs) plus a `metrics` table (the quantitative tuples), so many metrics attach to one coordinate without schema explosion.

`SET = (I, M, C, Ax, In, P)` as the coordinate; `Q` lives as the attached metric payload.

**Megatensor.** `MT = union over datasets d of SET_d`. Sparse by construction (only observed coordinates exist), append-only (new datasets contribute SETs, no refactor), columnar (Arrow/Parquet), and composably sliceable on any axis subset. Interoperability is implicit because all SETs share the axis system.

---

## 4. Canonical schema (the harmonization contract)

Every ingestion adapter, canon or PRIDE, emits rows in **one long-form observation schema**. This is the only boundary where per-source code lives. Everything downstream is generic.

### 4.1 Observation row (adapter output)

```text
# provenance / source
dataset_id            str        stable id you assign, e.g. "atlas4" / "oglcnacdb" / "PXD000000"
source_engine         enum       {fragpipe, maxquant, pd, diann, skyline, atlas, oglcnacdb}
source_file           str        relative path of the parsed file
file_checksum         str        sha256 of source_file

# identity (raw, pre-registry)
protein_id_raw        str
protein_id_type       enum       {uniprot_acc, uniprot_id, gene, refseq, protein_name}
isoform_raw           str|null
residue_pos_raw       int|null   position in the coordinate system of the source
residue_aa            str|null   S | T (others flagged)

# ptm
ptm_label             str        "O-GlcNAc"
ptm_unimod            str        "UniMod:43" (HexNAc, monoisotopic 203.079373). See §6 note.
loc_score             float|null
loc_method            enum       {ptmshepherd, ptmrs, mq_locprob, ascore, diann_ptm,
                                   manual_curation, atlas_unambiguous, atlas_ambiguous, none}
loc_is_ambiguous      bool

# condition
cond_cell_line        str|null
cond_tissue           str|null
cond_treatment        str|null
cond_timepoint        str|null
cond_fraction         str|null   sup | pel | whole | ...
cond_replicate        str|null
cond_replicate_type   enum|null  {biological, technical}

# acquisition
acq_ms_mode           enum|null  {DDA, DIA, PRM, unknown}
acq_msn_level         enum|null  {MS2, MS3, unknown}
acq_gradient_min      float|null
acq_collision         str|null   HCD | CID | ETD | EThcD | ...
acq_enrichment        str|null   lectin | antibody | chemoenzymatic(IsoTaG/etc) | metabolic | none

# instrument
inst_vendor           str|null
inst_model            str|null
inst_ms_cv            str|null   PSI-MS accession, e.g. MS:1002732 (Orbitrap Fusion Lumos)

# provenance (extended)
prov_pxd              str|null
prov_lab_pi           str|null
prov_country          str|null
prov_doi              str|null
prov_search_software  str|null
prov_search_version   str|null

# metric payload (one row per metric per event; melt wide quant here)
metric_name           enum       {intensity, abundance, ratio, log2fc, spectral_count,
                                   qvalue, score, quant_event}
metric_value          float
metric_norm_state     enum       {raw, normalized, log2, imputed, curated}
metric_unit           str|null
qc_flags              str|null   pipe-delimited
```

Wide quant columns (e.g. three intensity columns for three conditions) are **melted**: one observation row per (site, condition, metric). Never store wide.

### 4.2 Registry-resolved fields (filled by §7)

```text
protein_acc           str        canonical UniProt accession
isoform               str|null   normalized isoform token
residue_pos           int        position on UniProt canonical sequence
identity_id           uint64     hash(protein_acc, isoform, residue_pos, residue_aa)
ptm_id                uint64     hash(ptm_unimod)
condition_id          uint64     controlled-vocab id
acquisition_id        uint64
instrument_id         uint64
provenance_id         uint64
set_uid               uint64     hash("{identity_id}|{ptm_id}|{condition_id}|{acquisition_id}|{instrument_id}|{provenance_id}")
```

### 4.3 On-disk layout

```text
megatensor/
  registry/
    identity_dim.parquet        identity_id -> (protein_acc, isoform, residue_pos, residue_aa, protein_level_only:bool)
    ptm_dim.parquet
    condition_dim.parquet       controlled vocab + raw->id crosswalk
    acquisition_dim.parquet
    instrument_dim.parquet
    provenance_dim.parquet
  sets/  (hive partitioned by dataset_id)
    set_coordinates/            set_uid + all *_id
  metrics/ (hive partitioned by dataset_id)
    set_metrics/                set_uid, metric_name, metric_value, metric_norm_state, metric_unit, qc_flags
  enrichment/
    identity_enrichment.parquet identity_id -> disorder, domain, region, seq_window, ...
    set_enrichment.parquet      set_uid -> any context-conditioned features (optional)
  views/
    megatensor.sql              DuckDB view materializing the union join across the above
```

The Megatensor is a **view**, not a giant materialized table: a DuckDB join of `set_coordinates` + dims + `set_metrics` + enrichment. Append = write new partitions under `sets/` and `metrics/`. `set_uid` makes appends idempotent.

---

## 5. Data acquisition

### 5.1 Canon spine (build first, lowest risk)

**The O-GlcNAc Database (MCW).** https://www.oglcnac.mcw.edu/download/ ("Download O-GlcNAc reference datasets" + "Download O-GlcNAc literature metadata"). ~20,000 proteins, 48 species. Curated sites with evidence levels, species, references. Python helpers exist: https://github.com/Synthaze/utilsovs . `source_engine = oglcnacdb`, `loc_method = manual_curation`. License: non-profit free, **cite Wulff-Fuentes et al. 2021, PMID 33479245**.

**O-GlcNAcAtlas 4.0.** https://oglcnac.org/atlas/ . Dataset-I (unambiguous sites) and Dataset-II (ambiguous sites); >19,000 unambiguous, >11,000 ambiguous, >8,000 proteins, multi-species; ~5,400 quantification events on ~3,130 unique sites. Peptides are mapped to UniProtKB. `source_engine = atlas`. Map: Dataset-I -> `loc_is_ambiguous=false, loc_method=atlas_unambiguous`; Dataset-II -> `loc_is_ambiguous=true, loc_method=atlas_ambiguous`. The quantification events populate the Quant axis with `metric_name=quant_event, metric_norm_state=curated`. **Cite Ma et al. 2021 (Glycobiology) and Hou et al. 2025 (J Mol Biol, O-GlcNAcAtlas 4.0).**

**Canon modeling rules.**
- Entries with a localized site -> emit SETs.
- Protein-level-only entries (no residue) -> record in `identity_dim` with `protein_level_only=true`, do **not** emit a SET. A SET requires a site. This is deliberate; say so if asked.
- Canon fills Identity and Provenance hard, and (Atlas 4.0) some Quant. It has essentially no Condition / Acquisition / Instrument. That gap is the point: canon supplies the identity backbone, fresh PRIDE data supplies experimental context. Do not present the missing canon quant as a defect.

**Canon-vs-canon as a result.** The two reference DBs do not fully agree on sites once mapped to UniProt-canonical coordinates. Their overlap and disagreement is the first interoperability figure: the Megatensor reconciling the two reference sources nobody else reconciles. Compute it.

### 5.2 Fresh PRIDE (build second, parallelizable)

Tool: `pride-ingest` (https://github.com/filiprumenovski/pride-ingest). It lands PRIDE Archive **metadata** (projects, files) into raw -> bronze -> silver Parquet and exposes instruments, countries, softwares, PTM strings, keywords, and the file download locations. It is the discovery/audit surface, not a result-table fetcher.

**Scope rule:** ingest only **deposited result tables**. Do **not** re-search raw files; that is not an overnight operation.

**Step 1: snapshot metadata.**
```bash
cp .env.example .env
pip install -e .[dev]
pride-ingest ingest --mode both --output-root ./data --snapshot-date $(date +%F)
pride-ingest build-silver --output-root ./data
```
Use `--no-verify-ssl` if EBI trust-chain issues appear. A `--sample-size` dev run is fine for a dry pass; the full snapshot is needed for real discovery.

**Step 2: find glyco projects spread across the axes that differ.** Query bronze with DuckDB (PTM strings and softwares live in bronze structs, not silver):
```sql
WITH p AS (
  SELECT accession, title, countries, instruments, softwares, keywords, identifiedPTMStrings
  FROM read_parquet('data/bronze/projects/snapshot_date=*/**/*.parquet', hive_partitioning=true)
)
SELECT
  accession, title, countries,
  list_transform(instruments, x -> x.name) AS instruments,
  list_transform(softwares,   x -> x.name) AS softwares
FROM p
WHERE list_contains(list_transform(keywords, k -> lower(k)), 'o-glcnac')
   OR lower(title) LIKE '%glcnac%'
   OR EXISTS (SELECT 1 FROM UNNEST(identifiedPTMStrings) t(ptm)
              WHERE lower(ptm.name) LIKE '%glcnac%' OR lower(ptm.name) LIKE '%hexnac%')
ORDER BY accession;
```

**Step 3: choose 4-5 for maximum spread.** Selection criteria, in priority order:
- Distinct **search engines** across the picks (FragPipe, MaxQuant, Proteome Discoverer, DIA-NN, Skyline). The engine spread drives which adapters you build (§6). Build adapters only for the engines your picks actually use. Nothing more.
- Distinct **instruments** and **countries**. Heterogeneity is the variable you are claiming to defeat; the diversity is the figure.
- **At least one conditioned dataset** (>= 2 conditions/timepoints/treatments). Figure C (trajectory) must come from here. Canon cannot supply it.

**Step 4: resolve download URLs.** From silver `project_files`, `primary_public_location_value` holds the FTP/HTTPS path; filter to result categories and engine result-file patterns:
```sql
SELECT project_accession, file_name, primary_public_location_value
FROM read_parquet('data/silver/project_files/snapshot_date=*/*.parquet', hive_partitioning=true)
WHERE project_accession IN ('PXD......','PXD......')   -- the chosen picks
  AND (file_category_name IN ('RESULT','SEARCH') OR lower(file_name) LIKE '%result%');
```
Download those files only. PRIDE base for direct fetch: `https://ftp.pride.ebi.ac.uk/pride/data/archive/<YYYY>/<MM>/<PXD>/<file>`.

---

## 6. Ingestion adapters (TSV/CSV -> observation rows)

One adapter per search engine, **built only for the engines the picks use**. Each adapter knows four things: the site-identifying columns, the modification column, the quant columns to melt, and the localization-score column + its method tag.

| engine | canonical result file | site cols | mod column | quant cols (melt) | loc score / method |
|---|---|---|---|---|---|
| FragPipe (MSFragger + PTM-Shepherd / IonQuant) | `psm.tsv`, `combined_modified_peptide.tsv` | Protein, Protein ID (UniProt), Peptide, Assigned Modifications | `Assigned Modifications` / `Observed Modifications` (HexNAc 203.0794) | Intensity columns per run/experiment | PTM-Shepherd / MSFragger localization -> `ptmshepherd` |
| MaxQuant | `HexNAc (ST)Sites.txt` or `GlcNAc (ST)Sites.txt` (+ `evidence.txt`) | Proteins, Positions within proteins, Amino acid | site table is mod-specific | `Intensity ...`, `Ratio ...` columns | `Localization prob` -> `mq_locprob` |
| Proteome Discoverer | exported PSMs / PeptideGroups `.txt`/`.xlsx` | Master Protein Accessions, Positions in Master Proteins | Modifications | Abundance columns | ptmRS / `ptmRS: Best Site Probabilities` -> `ptmrs` |
| DIA-NN | `report.tsv` / `report.parquet` | Protein.Group, Protein.Ids, Modified.Sequence | within Modified.Sequence (UniMod:43) | Precursor/Fragment Quantity, MaxLFQ | DIA-NN PTM scoring -> `diann_ptm` |
| Skyline | exported `.csv` (Document Grid / pivot) | Protein, Peptide, modification annotations | modification annotation col | Area / Normalized Area | none/manual -> `none` |

Adapter responsibilities:
- Parse only the columns above; ignore everything else.
- Melt wide quant to one observation row per (site, condition, metric).
- Map condition/acquisition/instrument from the PRIDE metadata you already pulled for that PXD (join on `dataset_id`), not from guesswork inside the result file.
- Emit `protein_id_type` honestly per engine (FragPipe/DIA-NN usually UniProt; MaxQuant/PD vary; Skyline often protein name).
- Stamp `source_file`, `file_checksum`, `prov_search_software`, `prov_search_version`.

---

## 7. Registry (identity normalization)

Resolve every observation to a canonical biological identity. Pluggable resolver interface so the existing UniProt ETL can back it; fall back to UniProt REST.

1. Map `protein_id_raw` + `protein_id_type` -> canonical UniProt accession (`protein_acc`). Track isoforms explicitly.
2. Map `residue_pos_raw` (source coordinates) -> `residue_pos` on the UniProt canonical sequence. Where the source used a peptide-local or isoform coordinate, remap. Flag and quarantine any site whose `residue_aa` does not match the canonical sequence at `residue_pos` (`qc_flags += seq_mismatch`); do not silently keep mismatches.
3. Confirm `residue_aa in {S, T}` for O-GlcNAc; flag others.
4. Emit `identity_id = hash(protein_acc, isoform, residue_pos, residue_aa)`.

Resolver source: UniProt REST `https://rest.uniprot.org/uniprotkb/{acc}.json` (sequence + features), or the local ETL. Cache resolutions to Parquet; the registry is reused across all datasets.

---

## 8. Context encoding (controlled vocabularies)

For Condition, Acquisition, Instrument, Provenance: normalize raw strings into controlled-vocab tokens, then assign integer ids in dim tables. Keep a raw->token->id crosswalk in each dim so nothing is lost and provenance is auditable.

- **Condition:** normalize cell line (e.g. "HEK293T" / "HEK 293 T" -> `HEK293T`), treatment, timepoint, fraction.
- **Acquisition:** DDA/DIA, MSn level, collision, enrichment chemistry, gradient.
- **Instrument:** vendor + model + PSI-MS CV accession. Pull instrument straight from `pride-ingest` `project_instruments` for PRIDE sources.
- **Provenance:** PXD, lab PI, **country**, DOI, search software + version. Pull from `pride-ingest` silver/bronze. **These must land as real, queryable axis values, not a footnote.** The headline query (five labs, five instruments, one site) depends on it.

Null is a first-class value. Missing context (e.g. canon has no instrument) stays null and shows up honestly in axis completeness.

---

## 9. SET assembly and Megatensor append

1. Join observation rows to the dims to attach `*_id`.
2. Compute `set_uid`. Collapse rows sharing a `set_uid` into one SET coordinate; their metrics become the attached `set_metrics` rows.
3. Write `sets/set_coordinates` and `metrics/set_metrics` partitioned by `dataset_id`. Idempotent on `set_uid`.
4. The Megatensor view (`views/megatensor.sql`) joins coordinates + dims + metrics (+ enrichment).

**Refinement vs the model doc (state this deliberately):** the model doc puts "localization confidence" inside the PTM axis. Keep only the modification **class** in the identity coordinate (`ptm_id = hash(UniMod:43)`). Attach `loc_score` + `loc_method` as **payload** on the SET, not as part of the coordinate. Localization confidence is source-dependent and not a stable identity; folding it into the coordinate would shatter joinability across engines. This is the correct refinement and a good talking point.

---

## 10. Enrichment (every field, local-first)

Enrichment attaches as tensor-side payloads keyed on `identity_id` (sequence/structure/domain) or `set_uid` (context-conditioned). **Two clean enrichments fully wired beat six half-wired ones that rate-limit out at 2am.** Local-first. Only enrich identities/sites that are actually in the Megatensor.

| field | source | local? | scope | output columns |
|---|---|---|---|---|
| **Intrinsic disorder** | metapredict (Holehouse lab), `pip install metapredict` | yes, CPU | all proteins in MT | `disorder_score` (per-site), `disorder_region:bool` |
| **Domain / region** | UniProt features (local ETL or `rest.uniprot.org/.../{acc}.json`) | semi | all proteins in MT | `region_type` (domain/disordered/linker/repeat), `domain_name` |
| **Sequence window** | UniProt sequence | yes | all sites | `seq_window` (k-residue window centered on site, k=±7 default) |
| **Pathway / GO enrichment** | gseapy Enrichr, `pip install gseapy` | API, one call | one gene list (the MT proteome) | one enrichment table + one figure (`GO_Biological_Process`, `KEGG`) |
| **SASA / structure** | AlphaFold DB PDB + `freesasa` (`pip install freesasa`) | yes once PDB pulled | **visualized sites only** (Figure C site + neighbors) | `sasa`, `plddt` (from B-factor) |
| **Co-occurring PTMs / motif** (optional) | within-MT co-membership; k-mer window | yes | optional | `coptm_count`, `motif_window` |
| **Conservation** (optional, cut if tight) | phyloP/phastCons | API | optional | `phylop` |

Scope guards:
- AlphaFold model URL: `https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb`. Full-proteome SASA overnight is a trap. Scope SASA to the sites you actually plot.
- gseapy: exactly one gene list, one figure. Cheap, biological, done.
- Disorder via metapredict, not an external API. Do not hit IUPred at 2am.
- Domains from the existing UniProt ETL come effectively free; prefer it over live REST for volume.

---

## 11. What to produce

### 11.1 Figure A: axis completeness
Per source, percent non-null across the seven axes, plus SET / protein / site / metric counts. The story the figure tells: canon fills Identity + Provenance (+ some Quant from Atlas 4.0), PRIDE fills Quant + Condition + Acquisition + Instrument. The Megatensor is the sparse union of both. Stacked or small-multiple bars.

### 11.2 Figure B: cross-dataset interoperability (UpSet)
Site-membership UpSet across all sources after mapping to UniProt-canonical coordinates. Two layers:
- canon-vs-canon (O-GlcNAc DB vs Atlas Dataset-I vs Dataset-II): reconciliation of the reference sources.
- the 4-5 PRIDE datasets: heterogeneous sources (different instrument/country/engine) collapsing into one coordinate space.
This is the slide Charlie remembers. Heterogeneity is the headline.

### 11.3 Figure C: single-site multi-condition trajectory
Pick one biologically meaningful, well-localized O-GlcNAc site present in a conditioned PRIDE dataset. Plot its metric across conditions/timepoints, pulled live from one Megatensor slice. Proves §6.3 of the proposal.

### 11.4 Live query demo (DuckDB)
A short `queries.sql` run at the bench. Required headline query: same O-GlcNAc site, every source/lab/instrument that observed it, in one statement. Plus: per-axis slicing examples (per-protein, per-instrument, per-condition) to prove composable slicing live.

### 11.5 Tensor exports (no modeling)
- Dense `site x condition` matrix for a chosen metric -> Parquet + `.npy`.
- Enriched `site x feature` matrix (disorder, region one-hot, window-derived features) -> Parquet + `.npy`.
- Print shapes and write a short `EXPORTS.md` describing how a downstream model would consume them. Do not train anything.

### 11.6 Report + repo
`report.md`: abstract, ontology, methods, the four figures with captions, the completeness table, the honest-framing paragraph, citations. `README.md`: one-command reproduction. Keep the model-doc narrative; swap the "demonstration plan" language for realized results.

---

## 12. Run of show (autonomous loop + checkpoints)

Run autonomously between checkpoints. Stop and ping Filip only at these gates.

```
Phase 0  Canon spine
  - download O-GlcNAc DB + Atlas 4.0; adapters -> observation rows -> registry -> SETs
  - compute canon axis completeness + canon-vs-canon overlap
  >> CHECKPOINT 0: show canon completeness + overlap numbers. (canon-only is already a full deliverable)

Phase 1  PRIDE discovery
  - snapshot metadata via pride-ingest; run glyco discovery query; rank by engine/instrument/country spread
  >> CHECKPOINT 1: present the 4-5 candidate PXDs + their spread. FILIP CHOOSES THE FIVE. Do not download before approval.

Phase 2  PRIDE ingestion
  - download chosen result tables; build adapters ONLY for the engines those picks use; parse -> observation rows
  >> CHECKPOINT 2: parse QC (rows, sites, null rates, seq_mismatch count) per dataset before SET assembly.

Phase 3  Megatensor assembly
  - registry-resolve all; assemble SETs; append to MT; build the view
  >> CHECKPOINT 3: headline numbers (total SETs, sites, sources) + draft Figures A and B.

Phase 4  Enrichment + finalize
  - disorder + domain + window + GSEA; SASA scoped to Figure C site; Figure C trajectory; exports; report
  >> CHECKPOINT 4: enriched completeness + final figures + report for review.

Phase X  ML  ::  OUT OF SCOPE THIS RUN. Do not implement. Leave the export hooks (11.5) in place.
```

Sequencing rationale: canon first because everything hangs off it and a canon-only build is whole. PRIDE second, parallelizable, degrades gracefully if it slips. Enrichment third (garnish). A short night still lands something complete.

---

## 13. Guardrails (failure modes to refuse)

1. **Do not re-search raw MS files.** Ingest deposited result tables only.
2. **Adapter coverage follows the picks.** Build no parser for an engine you did not pull.
3. **Localization confidence is not comparable across sources.** Store `loc_score` + `loc_method` together; never unify into one number. FragPipe (PTM-Shepherd), MaxQuant (loc prob), PD (ptmRS), DIA-NN, and the curated DBs all mean different things by "confidence."
4. **Do not barf out every API.** Two clean enrichments over six flaky ones. Local-first. Cap external calls; only enrich identities in the MT.
5. **Never densify the full tensor.** Sparse Parquet + DuckDB view throughout. Dense only for the small explicit export matrices in 11.5.
6. **No ML this run.**
7. **Provenance axis carries instrument/country/software/PXD as real values.** That is the headline; do not bury it as metadata.
8. **Quarantine, do not silently drop.** Sites failing registry validation get `qc_flags`, not deletion. Report quarantine counts.

---

## 14. Repo layout and one-command run

```text
fehl-megatensor/
  pyproject.toml
  README.md
  justfile                 # or Makefile: `just all`
  src/megatensor/
    schema.py              # the §4 contract + dtypes
    ingest/
      base.py              # adapter interface -> observation rows
      fragpipe.py maxquant.py pd.py diann.py skyline.py   # only those needed
      canon_oglcnacdb.py canon_atlas.py
    registry.py            # §7
    context.py             # §8 controlled vocab + dims
    sets.py                # §9 SET assembly + append
    enrich/
      disorder.py domain.py window.py gsea.py sasa.py
    viz/
      completeness.py upset.py trajectory.py
    export.py              # §11.5
    report.py
  data/                    # pride-ingest output root + canon downloads (gitignored)
  megatensor/              # §4.3 on-disk layout (gitignored)
  views/megatensor.sql
  queries.sql              # §11.4
  report.md
```

```bash
just all   # phase0 -> phase4, stopping at checkpoints
# or granular:
just canon          # phase 0
just pride-discover # phase 1
just pride-ingest   # phase 2 (after approval)
just assemble       # phase 3
just enrich figures # phase 4
```

---

## 15. Data sources and citations

- **The O-GlcNAc Database (MCW).** https://www.oglcnac.mcw.edu/download/ . Cite: Wulff-Fuentes E, et al. The human O-GlcNAcome database and meta-analysis. *Sci Data* 8, 25 (2021). PMID 33479245. Tools: https://github.com/Synthaze/utilsovs . Non-profit use only; acknowledge in publications.
- **O-GlcNAcAtlas 4.0.** https://oglcnac.org/atlas/ . Cite: Ma J, et al. O-GlcNAcAtlas. *Glycobiology* 31(7):719-723 (2021), DOI 10.1093/glycob/cwab003; Hou C, et al. O-GlcNAcAtlas 4.0. *J Mol Biol* (2025), 169033.
- **PRIDE Archive.** https://www.ebi.ac.uk/pride/ ; v3 API under `/pride/ws/archive/v3/`. Accessed via `pride-ingest` (https://github.com/filiprumenovski/pride-ingest).
- **UniProt.** https://rest.uniprot.org/ . **AlphaFold DB.** https://alphafold.ebi.ac.uk/ .
- **metapredict.** https://github.com/idptools/metapredict . **gseapy.** https://github.com/zqfang/GSEApy . **freesasa.** https://freesasa.github.io/ .
- O-GlcNAc as HexNAc on Ser/Thr: **UniMod:43**, monoisotopic delta 203.079373. Match the engine's HexNAc/O-GlcNAc modification accordingly.
```
