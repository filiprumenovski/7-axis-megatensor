"""Engine-native score/quant columns -> long-form metric rows."""

from __future__ import annotations

import re
from typing import Any

MetricRow = tuple[str, float, str, str | None]  # name, value, norm_state, unit

FLOAT_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")


def _fval(raw: Any) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        v = float(raw)
        return v if v == v else None  # NaN check
    except (TypeError, ValueError):
        return None


def mq_metrics_for_row(row: dict, quant_col: str) -> list[MetricRow]:
    val = _fval(row.get(quant_col))
    if val is not None and val > 0:
        return [("intensity", val, "raw", "MQ_intensity")]
    return []


def pd_metrics_for_row(row: dict, cols: dict[str, str]) -> list[MetricRow]:
    out: list[MetricRow] = []
    intensity_col = cols.get("intensity")
    if intensity_col:
        val = _fval(row.get(intensity_col))
        if val is not None:
            out.append(("intensity", val, "raw", "PD_intensity"))
    qvals: list[float] = []
    for col_key, unit in (
        ("pep 1d", "PEP"),
        ("pep 2d", "PEP2D"),
        ("percolator q-value", "Percolator"),
    ):
        col = cols.get(col_key)
        if not col:
            continue
        val = _fval(row.get(col))
        if val is not None:
            qvals.append(val)
    if qvals:
        out.append(("qvalue", min(qvals), "raw", "PEP"))
    if not out:
        out.append(("spectral_count", 1.0, "raw", "PSM"))
    return out


def mztab_metrics_for_row(row: dict) -> list[MetricRow]:
    pep = _fval(row.get("search_engine_score[1]"))
    if pep is not None:
        return [("qvalue", pep, "raw", "mzTab_PEP"), ("spectral_count", 1.0, "raw", "PSM")]
    return [("spectral_count", 1.0, "raw", "PSM")]


def append_metric_rows(
    rows: list[dict],
    base: dict,
    metrics: list[MetricRow],
) -> None:
    if not metrics:
        metrics = [("spectral_count", 1.0, "curated", None)]
    for metric_name, metric_value, norm_state, unit in metrics:
        rows.append(
            {
                **base,
                "metric_name": metric_name,
                "metric_value": metric_value,
                "metric_norm_state": norm_state,
                "metric_unit": unit,
            }
        )
