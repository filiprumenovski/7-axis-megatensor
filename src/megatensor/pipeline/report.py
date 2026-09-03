"""Generate report.md from pipeline summaries."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from megatensor.paths import FIGURES, ROOT
from megatensor.store import CANON_STORE, PRIDE_STORE, UNION_STORE


def _fmt_n(v) -> str:
    if isinstance(v, (int, float)):
        return f"{int(v):,}"
    return str(v)


def _figure_c_label() -> str:
    traj_path = FIGURES / "figure_c_trajectory_candidate.csv"
    if not traj_path.is_file():
        return "—"
    traj = pl.read_csv(traj_path)
    if traj.is_empty():
        return "—"
    row = traj.row(0, named=True)
    return f"{row['protein_id_raw']}:{row['residue_pos_raw']}:{row['residue_aa']} — Light/Heavy SILAC in {row['dataset_id']}"


def _panel_figures() -> str:
    manifest = FIGURES / "panel_figures.json"
    if not manifest.is_file():
        return "_Run `megatensor figures` (requires `cmcrameri`)._"
    paths = json.loads(manifest.read_text())
    lines = []
    for key, exports in paths.items():
        pdf = Path(exports.get("pdf", "")).name if isinstance(exports, dict) else ""
        png = Path(exports.get("png", "")).name if isinstance(exports, dict) else Path(exports).name
        if pdf:
            lines.append(f"- **{key}:** `figures/{pdf}` (vector) · `figures/{png}` (preview)")
        else:
            lines.append(f"- **{key}:** `figures/{png}`")
    return "\n".join(lines)


def run_report() -> Path:
    canon = json.loads(CANON_STORE.summary_path.read_text()) if CANON_STORE.summary_path.is_file() else {}
    pride = json.loads(PRIDE_STORE.summary_path.read_text()) if PRIDE_STORE.summary_path.is_file() else {}
    union = (
        json.loads(UNION_STORE.summary_path.read_text()) if UNION_STORE.summary_path.is_file() else {}
    )
    enrich = (
        json.loads((UNION_STORE.enrichment / "enrichment_summary.json").read_text())
        if (UNION_STORE.enrichment / "enrichment_summary.json").is_file()
        else {}
    )

    fig_c = _figure_c_label()
    comp = enrich.get("completeness", {})
    gsea = enrich.get("gsea")
    sasa = enrich.get("figure_c_sasa")

    gsea_line = (
        f"Pathway enrichment (Enrichr): {gsea['terms']} terms from {gsea['genes']} genes → `{Path(gsea['path']).name}`"
        if gsea
        else "_GSEA skipped (install `gseapy` or network unavailable)._"
    )
    sasa_line = (
        f"Figure C SASA (AlphaFold + freesasa): **{sasa['sasa_total']} Å²** at {sasa['protein_acc']}:{sasa['residue_pos']} "
        f"({sasa.get('note', '')})"
        if sasa
        else "_SASA skipped (install `freesasa` or AlphaFold model unavailable)._"
    )

    human_shared = union.get("human_shared_sites", "—")
    human_pride_only = union.get("human_pride_only_sites", "—")

    text = f"""# 7-Axis Megatensor — Technical Report

## Panel headline

| | Canon | PRIDE |
|---|------:|------:|
| Observation rows | {_fmt_n(canon.get('observation_rows', '—'))} | {_fmt_n(pride.get('observation_rows', '—'))} |
| Unique sites | {_fmt_n(canon.get('unique_sites', '—'))} | {_fmt_n(pride.get('unique_sites', '—'))} |
| SETs | {_fmt_n(canon.get('unique_sets', '—'))} | {_fmt_n(pride.get('unique_sets', '—'))} |

**Cross-layer:** {_fmt_n(union.get('shared_sites', '—'))} shared sites · {_fmt_n(union.get('pride_only_sites', '—'))} PRIDE-only · {_fmt_n(union.get('canon_only_sites', '—'))} canon-only

**Human PRIDE** (rice PXD036527 excluded): {_fmt_n(human_shared)} shared · {_fmt_n(human_pride_only)} PRIDE-only novel vs canon

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

{_panel_figures()}

### Figure A — Axis completeness

Canon fills identity and reference provenance; PRIDE fills condition, acquisition, instrument, and engine-native quant metrics. CSV: `figures/axis_completeness_*.csv`.

### Figure B — Cross-layer interoperability

UpSet membership tables: `figures/upset_canon_membership.csv`, `figures/upset_pride_membership.csv`.
Human-filtered overlap: `figures/canon_vs_pride_overlap_human.json`.

### Figure C — Single-site trajectory

**{fig_c}**. Mean Light vs Heavy intensities in `figures/figure_c_trajectory.png`.
{sasa_line}

## Enrichment

| Feature | Coverage |
|---------|----------|
| Seq window | {comp.get('seq_window_pct', '—')}% |
| Domain/region | {comp.get('domain_pct', '—')}% |
| Disorder (metapredict) | {comp.get('disorder_pct', '—')}% |
| Gene symbol | {comp.get('gene_symbol_pct', '—')}% |

{gsea_line}

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
"""
    out = ROOT / "report.md"
    out.write_text(text)
    return out
