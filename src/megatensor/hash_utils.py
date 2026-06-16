"""Stable uint64 hashes for registry and SET coordinates."""

from __future__ import annotations

import xxhash


def hash_u64(*parts: str | int | None) -> int:
    """Hash string parts into a stable unsigned 64-bit integer."""
    h = xxhash.xxh64()
    for p in parts:
        h.update(b"|")
        if p is None:
            h.update(b"\x00")
        else:
            h.update(str(p).encode("utf-8"))
    return h.intdigest()


PTM_UNIMOD = "UniMod:43"
PTM_ID = hash_u64(PTM_UNIMOD)
