"""Intrinsic disorder via metapredict (optional dependency)."""

from __future__ import annotations

import structlog

log = structlog.get_logger()


def disorder_at_site(sequence: str, pos: int) -> float | None:
    """Per-residue disorder score at 1-based position, or None if unavailable."""
    try:
        import metapredict
    except ImportError:
        return None
    if not sequence or pos < 1 or pos > len(sequence):
        return None
    try:
        scores = metapredict.predict_disorder(sequence)
        return float(scores[pos - 1])
    except Exception as exc:
        log.debug("disorder_failed", error=str(exc))
        return None
