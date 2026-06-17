# Why the Megatensor Matters

*Defense / slide narrative — tied to pipeline outputs in `figures/megatensor_importance.json`.*

---

## The problem (one sentence)

O-GlcNAc proteomics is published as **flat, incompatible tables** that mix **what** was modified (site identity) with **how** it was measured (engine, tissue, SILAC, probe)—so cross-lab science and ML stay stuck in spreadsheet hell.

---

## What we built

A **sparse Site Event Tensor (SET)** megatensor: seven axes that **separate identity from context**, union canon libraries + PRIDE deposits **append-only**, query with DuckDB/Polars, export ML-ready tensors.

**Honest scope:** thin per-file adapters at ingest (~dozen lines each). The hard part we automate is **harmonization into the axis system**, not reading CSVs.

---

## Three payoffs (with numbers from this repo)

### 1. Structural interoperability

| Flat-file meta-analysis | Megatensor |
|-------------------------|------------|
| **66** pairwise PXD reconciliations for 12 studies | **12** one-time adapters |
| Manual site matching per comparison | One `site_key` join |
| Canon vs PRIDE = separate projects | **4,376** shared sites automatically |
| “Is this site replicated?” = ad hoc | **353** triangulated (canon + ≥2 PXDs) |

**Slide line:** *We didn’t merge spreadsheets—we merged coordinate systems.*

### 2. Context-preserving biology

Same `site_key` carries:

- **393** SILAC-quantified sites (Heavy/Light axis)
- **811** brain–liver site pairs (tissue axis, BAP1KO)
- **4,622** chemoproteomic sites (probe axis)
- **82** sites that are SILAC-quantified **and** triangulated

**Slide line:** *Tissue, SILAC, and chemo are slices—not separate databases.*

### 3. Evidence ranking + ML readiness

- Replication tiers, protein hubs (**HCFC1**: 18 triangulated sites)
- Concordance only meaningful after identity alignment (GlycoID *r*≈0.95 vs cross-lab *r*≈0.27)
- Exports: `site_x_condition.parquet`, `site_x_features.parquet`

**Slide line:** *The tensor turns deposits into ranked evidence, not row archives.*

---

## What you can say you **cannot** do easily without it

1. Rank **353** triangulated sites across canon + 3 PXDs without N² manual joins  
2. Plot brain vs liver for **811** sites and SILAC FC for the **same keys** in one pipeline  
3. Quantify “within-pipeline vs cross-lab” concordance at shared sites across China/US engines  
4. One DuckDB query: all observations for `P51610:579:T` across PXD035902, PXD039536, PXD058744  

See `queries/queries.sql`.

---

## What we are **not** claiming

- Global intensity calibration across engines (concordance is weak cross-lab—that’s a **finding**)
- Unified localization scores (ptmRS ≠ MaxQuant loc prob)
- Chemistry resolution (HexNAc vs O-GlcNAc collapsed to UniMod:43)
- Exhaustive PRIDE coverage

**The claim:** **representation** enables the right cross-study questions; quant harmonization is future work.

---

## Suggested defense arc (5 min)

1. **Problem** — flat files entangle identity + context (30 s)  
2. **Idea** — SET / 7-axis megatensor (45 s)  
3. **Proof** — Figure 0 impact + 4,376 shared / 353 triangulated (60 s)  
4. **Biology** — SILAC, brain/liver, OGT concordance thread (90 s)  
5. **Limitations** — honest, then ML exports + community resource (45 s)  

**Key figure:** `figures/analysis_megatensor_impact.pdf`

---

## Reproduce

```bash
just analyze && just analysis-figures && just biorxiv
```
