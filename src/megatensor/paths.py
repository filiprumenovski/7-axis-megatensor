"""Repository path constants."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
CANON = DATA / "canon"
PRIDE = DATA / "pride"
MT = ROOT / "megatensor"
REGISTRY = MT / "registry"
SETS = MT / "sets" / "set_coordinates"
METRICS = MT / "metrics" / "set_metrics"
ENRICHMENT = MT / "enrichment"
FIGURES = ROOT / "figures"
EXPORTS = ROOT / "exports"
VIEWS = ROOT / "views"
CACHE = ROOT / ".cache"
