"""Solvent-accessible surface area at a site via AlphaFold structure + freesasa."""

from __future__ import annotations

import tempfile
from pathlib import Path

import httpx
import structlog

log = structlog.get_logger()

AF_URL = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb"


def _download_alphafold_pdb(uniprot_acc: str, dest: Path) -> bool:
    url = AF_URL.format(acc=uniprot_acc.upper())
    try:
        r = httpx.get(url, timeout=60, follow_redirects=True)
        if r.status_code != 200 or not r.text.startswith("HEADER"):
            return False
        dest.write_text(r.text)
        return True
    except Exception as exc:
        log.warning("alphafold_download_failed", acc=uniprot_acc, error=str(exc))
        return False


def sasa_at_site(uniprot_acc: str, position: int) -> dict | None:
    """Per-residue SASA (Å²) at 1-based UniProt position from AlphaFold model."""
    try:
        import freesasa
    except ImportError:
        log.info("freesasa_skip", reason="not installed")
        return None

    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = Path(tmp) / f"AF-{uniprot_acc}.pdb"
        if not _download_alphafold_pdb(uniprot_acc, pdb_path):
            return None
        try:
            structure = freesasa.Structure(str(pdb_path))
            result = freesasa.calc(structure)
        except Exception as exc:
            log.warning("freesasa_failed", acc=uniprot_acc, error=str(exc))
            return None

        # Match residue number in PDB (may differ from UniProt for signal peptides).
        best: tuple[float, int] | None = None
        for i in range(structure.nAtoms()):
            res_seq = structure.residueNumber(i)
            chain = structure.chainLabel(i)
            sasa = result.atomArea(i)
            if res_seq == position:
                total = (best[0] if best else 0.0) + sasa
                best = (total, res_seq)
        if best is None:
            # Fallback: closest residue within ±3
            candidates: dict[int, float] = {}
            for i in range(structure.nAtoms()):
                res_seq = structure.residueNumber(i)
                if abs(res_seq - position) <= 3:
                    candidates[res_seq] = candidates.get(res_seq, 0.0) + result.atomArea(i)
            if not candidates:
                return None
            res_seq = min(candidates, key=lambda r: abs(r - position))
            return {
                "protein_acc": uniprot_acc,
                "residue_pos": position,
                "pdb_residue": res_seq,
                "sasa_total": round(candidates[res_seq], 2),
                "sasa_relative": None,
                "source": "alphafold_v4",
                "note": "nearest_pdb_residue",
            }

        return {
            "protein_acc": uniprot_acc,
            "residue_pos": position,
            "pdb_residue": best[1],
            "sasa_total": round(best[0], 2),
            "sasa_relative": round(best[0] / 200.0, 3),  # crude vs max ~200 Å²
            "source": "alphafold_v4",
            "note": "exact_match",
        }
