"""On-disk SET store layout — canon and PRIDE are separate tensorizations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from megatensor.paths import MT


@dataclass(frozen=True)
class TensorStore:
    """One self-contained Megatensor partition (registry + sets + metrics)."""

    root: Path
    layer: str  # canon | pride | union

    @property
    def registry(self) -> Path:
        return self.root / "registry"

    @property
    def sets(self) -> Path:
        return self.root / "sets" / "set_coordinates"

    @property
    def metrics(self) -> Path:
        return self.root / "metrics" / "set_metrics"

    @property
    def staging(self) -> Path:
        return self.root / "staging"

    @property
    def enrichment(self) -> Path:
        return self.root / "enrichment"

    @property
    def summary_path(self) -> Path:
        return self.root / f"{self.layer}_summary.json"


CANON_STORE = TensorStore(MT / "canon", "canon")
PRIDE_STORE = TensorStore(MT / "pride", "pride")
UNION_STORE = TensorStore(MT / "union", "union")
