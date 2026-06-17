"""Import-time figure preamble — FIGURES.md §1."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from cycler import cycler  # noqa: E402
from matplotlib import rcParams  # noqa: E402

STYLE_PATH = Path(__file__).resolve().parent / "style" / "house.mplstyle"

rcParams["pdf.fonttype"] = 42
rcParams["ps.fonttype"] = 42
rcParams["svg.fonttype"] = "none"
rcParams["svg.hashsalt"] = "figures"
rcParams["axes.prop_cycle"] = cycler(
    color=["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
)

# Okabe-Ito + panel accents
COLOR_MUTED = "#D8D8D8"
COLOR_MUTED_TEXT = "#666666"
COLOR_HIGHLIGHT = "#0072B2"
COLOR_ACCENT = "#E69F00"
COLOR_POS = "#009E73"
COLOR_NEG = "#D55E00"
COLOR_CANON = "#4C72B0"
COLOR_PRIDE = "#C44E52"
COLOR_PD = "#D55E00"
COLOR_MQ = "#0072B2"

# Clean blue sequential for 0–100% fill (avoid muddy batlow on binary-ish data)
FILL_CMAP_COLORS = ("#ffffff", "#e8f1fa", "#9ecae1", "#3182bd", "#08519c")
