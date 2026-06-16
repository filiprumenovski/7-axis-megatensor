"""CLI entrypoints."""

from __future__ import annotations

import subprocess
import typer

from megatensor.paths import ROOT

app = typer.Typer(help="FEHL Megatensor pipeline")


@app.command()
def canon() -> None:
    """Phase 0: ingest canon bulk CSVs -> SETs."""
    from megatensor.pipeline.canon import run_canon

    run_canon()


@app.command("pride-discover")
def pride_discover(
    output_root: str = typer.Option("data/pride", help="pride-ingest output root"),
    sample_size: int | None = typer.Option(None, help="dev sample size"),
) -> None:
    """Phase 1: PRIDE metadata snapshot + glyco discovery query."""
    cmd = [
        "pride-ingest",
        "ingest",
        "--mode",
        "both",
        "--output-root",
        output_root,
        "--snapshot-date",
        subprocess.check_output(["date", "+%F"], text=True).strip(),
    ]
    if sample_size:
        cmd.extend(["--sample-size", str(sample_size)])
    subprocess.run(cmd, check=True, cwd=ROOT)
    subprocess.run(["pride-ingest", "build-silver", "--output-root", output_root], check=True, cwd=ROOT)
    typer.echo("Run glyco discovery with queries/pride_glyco_discovery.sql")


@app.command("pride-ingest")
def pride_ingest() -> None:
    typer.echo("Phase 2: awaiting Filip's PXD picks (checkpoint 1).")


@app.command()
def assemble() -> None:
    typer.echo("Phase 3: assemble Megatensor view (after PRIDE parse).")


@app.command()
def enrich() -> None:
    typer.echo("Phase 4: enrichment (not yet wired).")


@app.command()
def figures() -> None:
    typer.echo("Figures: run after enrich.")


@app.command()
def export() -> None:
    typer.echo("Tensor exports: run after enrich.")


@app.command()
def report() -> None:
    typer.echo("Report: run at end.")


if __name__ == "__main__":
    app()
