"""Repository path constants."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
CANON = DATA / "canon"
PRIDE = DATA / "pride"
PRIDE_DOWNLOADS = PRIDE / "downloads"
MT = ROOT / "megatensor"
FIGURES = ROOT / "figures"
EXPORTS = ROOT / "exports"
VIEWS = ROOT / "views"
CACHE = ROOT / ".cache"

# Legacy aliases — prefer megatensor.store.CANON_STORE / PRIDE_STORE
REGISTRY = MT / "canon" / "registry"
SETS = MT / "canon" / "sets" / "set_coordinates"
METRICS = MT / "canon" / "metrics" / "set_metrics"
ENRICHMENT = MT / "canon" / "enrichment"
