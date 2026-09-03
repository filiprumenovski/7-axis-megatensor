"""CLI entrypoints."""

from __future__ import annotations

import typer

app = typer.Typer(help="7-Axis Megatensor pipeline")


@app.command()
def canon() -> None:
    """Phase 0: canon reference libraries -> isolated canon tensor."""
    from megatensor.pipeline.canon import run_canon

    run_canon()


@app.command("pride-discover")
def pride_discover(
    output_root: str = typer.Option("data/pride", help="pride-ingest output root"),
    sample_size: int | None = typer.Option(None, help="Limit records for live ingest dev runs"),
    backend: str = typer.Option("api", help="pride-ingest backend when --live-ingest"),
    live_ingest: bool = typer.Option(False, help="Pull fresh PRIDE metadata instead of local snapshot"),
    snapshot_date: str | None = typer.Option(None, help="Override snapshot_date partition"),
) -> None:
    """Phase 1: glyco discovery from local snapshot (default) or live pride-ingest."""
    from pathlib import Path

    from megatensor.pipeline.pride_discover import run_pride_discover

    run_pride_discover(
        Path(output_root),
        sample_size=sample_size,
        backend=backend,
        live_ingest=live_ingest,
        snapshot_date=snapshot_date,
    )


@app.command("pride-download")
def pride_download(
    dest: str = typer.Option("data/pride/downloads", help="Download destination root"),
    accession: list[str] = typer.Option(None, "--accession", "-a", help="Limit to PXD(s)"),
    dry_run: bool = typer.Option(False, help="Print manifest only"),
) -> None:
    """Phase 2: download curated PRIDE result tables via Aspera."""
    from pathlib import Path

    from megatensor.ingest.pride_download import download_picks

    summary = download_picks(Path(dest), accessions=accession or None, dry_run=dry_run)
    if not dry_run:
        typer.echo(summary.group_by("status").len())


@app.command("pride-ingest")
def pride_ingest_tables(
    dest: str = typer.Option("data/pride/downloads", help="Download root with parsed files"),
) -> None:
    """Phase 2b: PRIDE adapters -> observation rows (no tensorization)."""
    from pathlib import Path

    from megatensor.pipeline.pride_ingest import run_pride_ingest

    run_pride_ingest(Path(dest))


@app.command("pride-tensorize")
def pride_tensorize(
    dest: str = typer.Option("data/pride/downloads", help="Download root"),
    reparse: bool = typer.Option(False, help="Re-parse files before tensorizing"),
) -> None:
    """Phase 3: PRIDE observations -> isolated experimental tensor."""
    from pathlib import Path

    from megatensor.pipeline.pride_tensorize import run_pride_tensorize

    run_pride_tensorize(Path(dest), reparse=reparse)


@app.command()
def union() -> None:
    """Cross-layer site index (canon vs PRIDE interoperability)."""
    from megatensor.pipeline.union import run_union

    run_union()


@app.command()
def assemble() -> None:
    """Deprecated alias for pride-tensorize."""
    typer.echo("assemble is deprecated — use: megatensor pride-tensorize")
    from megatensor.pipeline.pride_tensorize import run_pride_tensorize

    run_pride_tensorize()


@app.command()
def enrich() -> None:
    """Phase 4a: UniProt window/domain enrichment (+ optional disorder/GSEA)."""
    from megatensor.pipeline.enrich import run_enrich

    run_enrich()


@app.command()
def figures() -> None:
    """Figure A/B/C datasets — separate canon vs PRIDE narratives."""
    from megatensor.pipeline.figures import run_figures

    run_figures()


@app.command()
def export() -> None:
    """Phase 4b: ML-ready site x condition / site x feature matrices."""
    from megatensor.export import run_export

    run_export()


@app.command()
def analyze() -> None:
    """Cross-layer analyses → megatensor/analysis/*.parquet."""
    from megatensor.pipeline.analyze import run_analyze

    run_analyze()


@app.command()
def analysis_figures() -> None:
    """Render analysis-backed figures (requires analyze)."""
    from megatensor.viz.analysis_plots import render_analysis_figures

    paths = render_analysis_figures()
    typer.echo(f"Wrote {len(paths)} analysis figures to figures/")


@app.command()
def biorxiv() -> None:
    """Write biorxiv.md preprint draft from analysis summary."""
    from megatensor.pipeline.biorxiv import run_biorxiv

    path = run_biorxiv()
    typer.echo(f"Wrote {path}")


@app.command()
def publish() -> None:
    """Full bioRxiv path: analyze → analysis figures → enrich → export → biorxiv."""
    from megatensor.pipeline.analyze import run_analyze
    from megatensor.pipeline.biorxiv import run_biorxiv
    from megatensor.pipeline.enrich import run_enrich
    from megatensor.pipeline.supplementary import write_supplementary
    from megatensor.export import run_export
    from megatensor.viz.analysis_plots import render_analysis_figures

    run_analyze(skip_gsea=True)
    run_enrich()
    from megatensor.pipeline.analyze import run_gsea_step

    run_gsea_step()
    render_analysis_figures()
    run_export()
    write_supplementary()
    path = run_biorxiv()
    typer.echo(f"Publish bundle ready — see {path} and figures/analysis_*.pdf")


@app.command()
def report() -> None:
    """Write report.md from pipeline summaries."""
    from megatensor.pipeline.report import run_report

    path = run_report()
    typer.echo(f"Wrote {path}")


@app.command()
def panel() -> None:
    """Panel polish: union stats, figures (PNG), enrichment garnish, report."""
    from megatensor.pipeline.enrich import run_enrich
    from megatensor.pipeline.figures import run_figures
    from megatensor.pipeline.report import run_report
    from megatensor.pipeline.union import run_union

    run_union()
    run_figures()
    run_enrich()
    path = run_report()
    typer.echo(f"Panel ready — see figures/*.png and {path}")


if __name__ == "__main__":
    app()
