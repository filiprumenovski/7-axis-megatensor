"""Phase 1: PRIDE glyco discovery from local snapshot or live ingest."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import polars as pl
import structlog

from megatensor.paths import FIGURES, PRIDE, ROOT

log = structlog.get_logger()

SNAPSHOT_TARBALL = ROOT / "pride_snapshot_parquets_2026-04-07.tar.gz"
GLYCO_SQL = (ROOT / "queries" / "pride_glyco_discovery.sql").read_text()
CATALOG_SQL = (ROOT / "queries" / "pride_glyco_catalog.sql").read_text()

# Curated PXDs maximizing engine + instrument + country spread (checkpoint 1 draft).
RECOMMENDED_PICKS: dict[str, str] = {
    "PXD035902": "PD — US — BAP1KO multi-tissue PSM tables",
    "PXD039536": "MaxQuant — China — Light/Heavy site tables (Figure C)",
    "PXD058744": "MaxQuant — US — Fusion+Exploris — drug resistance sites",
    "PXD064117": "mzTab — US — FOXK2/OGT ferroptosis dual conditions",
    "PXD064782": "mzTab — China — LATS1 tumor suppression O-GlcNAcylation",
    "PXD033062": "MaxQuant zip — US — GlycoID serum 30min ST sites (proximity)",
    "PXD042838": "PD xlsx — US — palbociclib resistance site table",
    "PXD033043": "MaxQuant zip — US — GlycoID insulin 30min ST sites",
    "PXD036527": "MaxQuant — China — rice OGT site tables",
    "PXD063995": "mzTab — bioorthogonal probe benchmark (4 condition tables)",
    "PXD033026": "MaxQuant zip — US — GlycoID serum cytosolic ST sites",
    "PXD014785": "MaxQuant — China — HexNAcSTSites.txt",
}


def find_snapshot_date(pride_root: Path) -> str | None:
    """Pick the snapshot partition with the most bronze project parquet parts."""
    projects = pride_root / "bronze" / "projects"
    if not projects.is_dir():
        return None
    best: tuple[int, str] | None = None
    for part in projects.glob("snapshot_date=*"):
        if not part.is_dir():
            continue
        snap = part.name.replace("snapshot_date=", "")
        n = len(list(part.glob("*.parquet")))
        if best is None or n > best[0]:
            best = (n, snap)
    return best[1] if best else None


def ensure_pride_data(pride_root: Path) -> str:
    snap = find_snapshot_date(pride_root)
    if snap:
        log.info("pride_snapshot_local", root=str(pride_root), snapshot_date=snap)
        return snap

    if SNAPSHOT_TARBALL.is_file():
        script = ROOT / "scripts" / "unpack_pride_snapshot.sh"
        subprocess.run(["bash", str(script), str(SNAPSHOT_TARBALL)], check=True, cwd=ROOT)
        snap = find_snapshot_date(pride_root)
        if snap:
            return snap

    raise FileNotFoundError(
        f"No PRIDE snapshot under {pride_root}. "
        f"Drop {SNAPSHOT_TARBALL.name} in repo root or run pride-ingest ingest."
    )


def _require_pride_ingest() -> str:
    candidates: list[Path] = []
    if venv := os.environ.get("VIRTUAL_ENV"):
        candidates.append(Path(venv) / "bin" / "pride-ingest")
    candidates.append(Path(sys.argv[0]).resolve().parent / "pride-ingest")
    for path in candidates:
        if path.is_file():
            return str(path)
    if cmd := shutil.which("pride-ingest"):
        return cmd
    raise RuntimeError(
        "pride-ingest not installed. pip install -e '.[pride]' or use local snapshot only."
    )


def _sql_for_root(sql: str, pride_root: Path, snapshot_date: str) -> str:
    out = sql.replace("data/pride/", f"{pride_root.as_posix()}/")
    out = re.sub(
        r"snapshot_date=\*/\*\*/\*\.parquet",
        f"snapshot_date={snapshot_date}/*.parquet",
        out,
    )
    out = re.sub(
        r"snapshot_date=\*/\*\.parquet",
        f"snapshot_date={snapshot_date}/*.parquet",
        out,
    )
    return out


def _engine_summary(catalog: pl.DataFrame) -> pl.DataFrame:
    if catalog.is_empty() or "engine_kind" not in catalog.columns:
        return pl.DataFrame()
    return (
        catalog.group_by("engine_kind")
        .agg(
            pl.col("accession").n_unique().alias("projects"),
            pl.col("result_files").sum().alias("result_files"),
        )
        .sort("projects", descending=True)
    )


def _write_picks(catalog: pl.DataFrame) -> pl.DataFrame:
    if catalog.is_empty():
        return catalog
    picks = catalog.filter(pl.col("accession").is_in(list(RECOMMENDED_PICKS.keys())))
    if picks.is_empty():
        return picks
    picks = picks.with_columns(
        pl.col("accession").replace_strict(RECOMMENDED_PICKS).alias("pick_note")
    )
    _write_flat_csv(picks, FIGURES / "pride_glyco_picks.csv")
    summary = {
        "pick_count": picks["accession"].n_unique(),
        "engines": picks["engine_kind"].unique().sort().to_list(),
        "astral": int(picks.filter(pl.col("has_astral")).height > 0),
        "fragpipe_in_archive": int(
            catalog.filter(pl.col("engine_kind").is_in(["fragpipe", "fragpipe_zip"])).height
        ),
    }
    (FIGURES / "pride_glyco_picks_summary.json").write_text(json.dumps(summary, indent=2))
    return picks


def run_pride_discover(
    output_root: Path | None = None,
    *,
    sample_size: int | None = None,
    snapshot_date: str | None = None,
    backend: str = "api",
    live_ingest: bool = False,
) -> pl.DataFrame:
    out = output_root or PRIDE
    out.mkdir(parents=True, exist_ok=True)

    if live_ingest:
        snap = snapshot_date or date.today().isoformat()
        pride_ingest = _require_pride_ingest()
        ingest_cmd = [
            pride_ingest,
            "ingest",
            "--mode",
            "both",
            "--output-root",
            str(out),
            "--snapshot-date",
            snap,
            "--backend",
            backend,
        ]
        if sample_size is not None:
            ingest_cmd.extend(["--sample-size", str(sample_size)])
        subprocess.run(ingest_cmd, check=True, cwd=ROOT)
        subprocess.run(
            [pride_ingest, "build-silver", "--output-root", str(out), "--snapshot-date", snap],
            check=True,
            cwd=ROOT,
        )
    else:
        snap = snapshot_date or ensure_pride_data(out)

    FIGURES.mkdir(parents=True, exist_ok=True)

    hits = _run_duckdb(_sql_for_root(GLYCO_SQL, out, snap))
    if hits.height:
        _write_flat_csv(hits, FIGURES / "pride_glyco_candidates.csv")
        log.info("glyco_metadata_candidates", count=hits.height)

    catalog = _run_duckdb(_sql_for_root(CATALOG_SQL, out, snap))
    if catalog.height:
        _write_flat_csv(catalog, FIGURES / "pride_glyco_catalog.csv")
        log.info("glyco_file_catalog", rows=catalog.height, projects=catalog["accession"].n_unique())
        summary = _engine_summary(catalog)
        if summary.height:
            _write_flat_csv(summary, FIGURES / "pride_engine_summary.csv")
            print(summary)

    picks = _write_picks(catalog)
    if picks.height:
        log.info("glyco_picks", count=picks["accession"].n_unique())
        print(picks.select("accession", "engine_kind", "title", "pick_note").unique())

    return catalog if catalog.height else hits


def _run_duckdb(sql: str) -> pl.DataFrame:
    try:
        import duckdb

        return duckdb.connect().execute(sql).pl()
    except ImportError:
        pass
    duckdb_cli = shutil.which("duckdb")
    if not duckdb_cli:
        return pl.DataFrame()
    result = subprocess.run(
        [duckdb_cli, "-csv", "-c", sql],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        log.warning("duckdb_query_failed", stderr=result.stderr.strip())
        return pl.DataFrame()
    import io

    return pl.read_csv(io.StringIO(result.stdout))


def _write_flat_csv(df: pl.DataFrame, path: Path) -> None:
    flat = df
    for col, dtype in df.schema.items():
        if dtype == pl.List(pl.Utf8) or dtype.base_type() == pl.List:
            flat = flat.with_columns(pl.col(col).list.join("|").alias(col))
        elif dtype.base_type() == pl.Struct:
            flat = flat.with_columns(pl.col(col).cast(pl.Utf8))
        elif dtype == pl.Boolean:
            flat = flat.with_columns(pl.col(col).cast(pl.Utf8))
    flat.write_csv(path)
