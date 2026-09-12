"""Shared helpers used across the pipeline modules.

Small, dependency-light utilities that several modules had each re-implemented: the repo
path constants, physical constants, the config-file loader, manifest/result readers, an
XYZ writer, and a tolerant float coercion. Import these instead of copy-pasting them so the
paths and conventions stay in one place.

Kept intentionally free of heavy imports (no rdkit/ase/pyscf at module load) so importing
`redox.common` is cheap and cannot introduce import cycles.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
from types import ModuleType

# --- Repo layout (single source of truth for these paths) ---
ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "config"
LIBRARY = ROOT / "library"
CALCS = ROOT / "calcs"
UMA = CALCS / "uma"
DFT = CALCS / "dft"
RESULTS = ROOT / "results"

# --- Physical constants ---
HARTREE_EV = 27.211386245988
EV_KJ = 96.485            # eV -> kJ/mol
EV_MEV = 1000.0
KT_EV = 0.0256926         # k_B * T at 298.15 K, in eV


def load_config(name: str) -> ModuleType:
    """Import and return a `config/<name>.py` module by path.

    Config files are loaded by path (not as installed packages) so the pipeline works from a
    plain checkout. Callers that want a single attribute do `getattr(load_config(name), attr)`.
    """
    spec = importlib.util.spec_from_file_location(name, CONFIG_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_manifest() -> list[dict]:
    """Rows of `library/manifest.csv` as a list of dicts (build.py's output)."""
    with (LIBRARY / "manifest.csv").open() as f:
        return list(csv.DictReader(f))


def read_result(gid: str, state: str, root: Path = DFT) -> dict | None:
    """Parsed `<root>/<gid>/<state>/result.json`, or None if it does not exist.

    Defaults to the DFT calc tree; pass `root=UMA` for UMA results.
    """
    p = root / gid / state / "result.json"
    return json.loads(p.read_text()) if p.exists() else None


def write_xyz(atoms, path: Path, comment: str = "") -> None:
    """Write an ASE Atoms object to a plain .xyz file (creating parent dirs)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), comment]
    for s, p in zip(atoms.get_chemical_symbols(), atoms.positions):
        lines.append(f"{s} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    path.write_text("\n".join(lines) + "\n")


def to_float(x):
    """Best-effort float; None on empty/None/non-numeric input."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return None
