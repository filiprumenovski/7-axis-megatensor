"""Download PRIDE result tables via Aspera (reuses ci-cd_nextflow ascp + keys)."""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

import polars as pl
import structlog

from megatensor.paths import PRIDE, ROOT
from megatensor.pipeline.pride_discover import RECOMMENDED_PICKS, ensure_pride_data, find_snapshot_date

log = structlog.get_logger()

DEFAULT_ASPERA_ROOT = Path("/run/media/filip/Data/ci-cd_nextflow/ci-cd")
MANIFEST_SQL = (ROOT / "queries" / "pride_download_manifest.sql").read_text()


def aspera_root() -> Path:
    return Path(os.environ.get("MEGATENSOR_ASPERA_ROOT", DEFAULT_ASPERA_ROOT))


def ascp_binary() -> Path:
    override = os.environ.get("FRAGPIPE_ASCP_BIN", "").strip()
    if override:
        return Path(override)
    candidate = aspera_root() / "tools" / "vendor" / "aspera-cli" / "bin" / "ascp"
    if candidate.is_file():
        return candidate
    launcher = aspera_root() / "tools" / "aspera" / "fragpipe_ascp.py"
    if launcher.is_file():
        return launcher
    raise FileNotFoundError(
        f"No ascp found under {aspera_root()}. Set MEGATENSOR_ASPERA_ROOT or FRAGPIPE_ASCP_BIN."
    )


def ssh_key_path() -> Path:
    override = os.environ.get("FRAGPIPE_ASPERA_SSH_KEY", "").strip()
    if override:
        return Path(override)
    etc = aspera_root() / "tools" / "aspera" / "etc"
    for name in ("aspera_tokenauth_id_rsa", "asperaweb_id_dsa.openssh"):
        key = etc / name
        if key.is_file():
            return key
    raise FileNotFoundError(f"No PRIDE Aspera SSH key in {etc}")


def to_aspera_spec(location: str) -> str:
    """Normalize PRIDE FTP/aspera location to prd_ascp@fasp.ebi.ac.uk:path."""
    location = location.strip()
    if location.startswith("aspera://"):
        location = location[len("aspera://") :]
    if location.startswith("ftp://ftp.pride.ebi.ac.uk/"):
        return f"prd_ascp@fasp.ebi.ac.uk:{location[len('ftp://ftp.pride.ebi.ac.uk/') :]}"
    if location.startswith("ftp://"):
        # generic ftp -> fasp host
        return re.sub(r"^ftp://[^/]+/", "prd_ascp@fasp.ebi.ac.uk:", location)
    if "@" in location and ":" in location.split("@", 1)[1]:
        return location
    raise ValueError(f"Unrecognized PRIDE location: {location}")


def _sql_paths(pride_root: Path, snap: str) -> str:
    sql = MANIFEST_SQL.replace("data/pride/", f"{pride_root.as_posix()}/")
    return re.sub(r"snapshot_date=\*/\*\.parquet", f"snapshot_date={snap}/*.parquet", sql)


def build_manifest(pride_root: Path | None = None) -> pl.DataFrame:
    root = pride_root or PRIDE
    snap = find_snapshot_date(root)
    if not snap:
        snap = ensure_pride_data(root)
    import duckdb

    return duckdb.connect().execute(_sql_paths(root, snap)).pl()


def download_file(
    location: str,
    dest: Path,
    *,
    ascp: Path | None = None,
    key: Path | None = None,
    bandwidth: str | None = None,
) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        log.info("skip_exists", dest=str(dest))
        return

    partial = dest.with_suffix(dest.suffix + ".partial")
    partial.unlink(missing_ok=True)

    ascp_path = ascp or ascp_binary()
    key_path = key or ssh_key_path()
    bw = bandwidth or os.environ.get("FRAGPIPE_ASPERA_BANDWIDTH", "200m")
    port = os.environ.get("FRAGPIPE_ASPERA_PORT", "33001")
    spec = to_aspera_spec(location)

    if ascp_path.name == "fragpipe_ascp.py":
        cmd = ["python3", str(ascp_path)]
    else:
        cmd = [str(ascp_path)]

    cmd.extend(
        [
            "-QT",
            "-k1",
            "-P",
            port,
            "-l",
            bw,
            "-i",
            str(key_path),
            spec,
            str(partial),
        ]
    )

    env = os.environ.copy()
    lib = aspera_root() / "tools" / "vendor" / "aspera-cli" / "install" / "lib"
    if lib.is_dir():
        env["LD_LIBRARY_PATH"] = f"{lib}:{env.get('LD_LIBRARY_PATH', '')}".rstrip(":")

    log.info("ascp_download", spec=spec, dest=str(dest))
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, check=False)
    if proc.returncode != 0:
        partial.unlink(missing_ok=True)
        raise RuntimeError(
            f"ascp failed ({proc.returncode}): {(proc.stderr or proc.stdout)[-2000:]}"
        )
    if not partial.is_file() or partial.stat().st_size == 0:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"ascp produced empty file for {spec}")

    partial.rename(dest)


def download_picks(
    dest_root: Path | None = None,
    *,
    accessions: list[str] | None = None,
    dry_run: bool = False,
) -> pl.DataFrame:
    """Download curated PXD result tables to data/pride/downloads/{PXD}/."""
    out_root = dest_root or (PRIDE / "downloads")
    manifest = build_manifest()
    if accessions:
        manifest = manifest.filter(pl.col("project_accession").is_in(accessions))

    log.info("download_manifest", files=manifest.height, projects=manifest["project_accession"].n_unique())
    if dry_run:
        print(manifest)
        return manifest

    results: list[dict] = []
    for row in manifest.iter_rows(named=True):
        pxd = row["project_accession"]
        fname = row["file_name"]
        dest = out_root / pxd / fname
        try:
            download_file(row["primary_public_location_value"], dest)
            results.append({"accession": pxd, "file_name": fname, "status": "ok", "path": str(dest)})
        except Exception as exc:
            log.error("download_failed", pxd=pxd, file=fname, error=str(exc))
            results.append({"accession": pxd, "file_name": fname, "status": "error", "error": str(exc)})

    summary = pl.DataFrame(results)
    manifest_path = out_root / "download_manifest.parquet"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_parquet(manifest_path)
    summary.write_csv(out_root / "download_results.csv")
    ok = summary.filter(pl.col("status") == "ok").height
    log.info("download_complete", ok=ok, failed=summary.height - ok)
    return summary
