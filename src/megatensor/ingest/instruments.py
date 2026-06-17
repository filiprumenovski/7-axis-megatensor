"""PSI-MS CV accession lookup for common instruments in our picks."""

from __future__ import annotations

# subset of PSI-MS; extend as needed
INSTRUMENT_CV: dict[str, str] = {
    "LTQ Orbitrap Elite": "MS:1000031",
    "LTQ Orbitrap Velos": "MS:1001742",
    "Orbitrap Fusion": "MS:1002416",
    "Orbitrap Fusion Lumos": "MS:1002732",
    "Orbitrap Exploris 480": "MS:1003095",
    "Orbitrap Astral": "MS:1004092",
    "Q Exactive HF": "MS:1002523",
    "Dionex UltiMate 3000 RSLCnano": "MS:1000526",
}


def psi_ms_cv(model: str | None) -> str | None:
    if not model:
        return None
    if model in INSTRUMENT_CV:
        return INSTRUMENT_CV[model]
    for key, cv in INSTRUMENT_CV.items():
        if key in model:
            return cv
    return None
