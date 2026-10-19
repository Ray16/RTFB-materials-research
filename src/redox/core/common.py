"""Shared helpers used across the pipeline modules.

Small, dependency-light utilities that several modules had each re-implemented: the repo
path constants, physical constants, the config-file loader, manifest/result readers, an
XYZ writer, and a tolerant float coercion. Import these instead of copy-pasting them so the
paths and conventions stay in one place.

Kept intentionally free of heavy imports (no rdkit/ase/pyscf at module load) so importing
`redox.core.common` is cheap and cannot introduce import cycles.
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path
from types import ModuleType

# --- Repo layout (single source of truth for these paths) ---
ROOT = Path(__file__).resolve().parents[3]
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


# energy fields that come ONLY from the active-protocol record (redox.core.protocol), never from
# result.json — see read_result.
ENERGY_KEYS = ("e_smd_Ha", "e_smd_eV", "converged_smd", "e_gas_Ha", "e_gas_eV",
               "converged_gas", "gas_homo_eV", "gas_lumo_eV", "gas_s_squared",
               "gas_spin_contam", "anion_unbound", "dG_solv_eV",
               "smd_homo_eV", "smd_lumo_eV", "smd_s_squared", "smd_spin_contam")


def read_result(gid: str, state: str, root: Path = DFT, raw: bool = False) -> dict | None:
    """Parsed `<root>/<gid>/<state>/result.json`, or None if it does not exist.

    For the DFT tree the ENERGY fields (ENERGY_KEYS) are replaced by the converged
    ACTIVE-protocol record (redox.core.protocol: uniform def2-TZVPD, hashed on geometry + level).
    If that record is missing or unconverged the energy fields are None and
    `sp_status == "missing"` — callers must treat the state as INCOMPLETE rather than use an
    energy computed at another level. Geometry/thermal/charge/spin fields come from
    result.json. `raw=True` returns the file untouched (provenance/audit only); pass
    `root=UMA` for UMA results.
    """
    p = root / gid / state / "result.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    if raw or root != DFT:
        return r
    from redox.core.protocol import load_record
    rec = load_record(p.parent, int(r["charge"]), int(r["mult"]))
    for k in ENERGY_KEYS:
        r[f"raw_{k}"] = r.get(k)
        r[k] = rec.get(k) if rec else None
    r["sp_status"] = "ok" if rec else "missing"
    r["sp_protocol_hash"] = rec.get("protocol_hash") if rec else None
    return apply_thermal_policy(r, p.parent)


def apply_thermal_policy(r: dict, state_dir: Path) -> dict:
    """Apply qualified native thermochemistry and an explicitly selected spin ensemble."""
    from redox.validation.spin_ensemble import apply_selected_ensemble
    return apply_selected_ensemble(_apply_single_surface_thermal_policy(r, state_dir), state_dir)


def _apply_single_surface_thermal_policy(r: dict, state_dir: Path) -> dict:
    """Overlay the explicitly selected thermal definition without a legacy fallback."""
    r = dict(r)
    gid, state = Path(state_dir).parent.name, Path(state_dir).name
    row = next((row for row in read_manifest() if row['id']==gid and row['state']==state), None)
    r['geometry_spin_qc'] = 'ok'
    if row is not None:
        from redox.qm.dft import _uma_mult
        try:
            desired = _uma_mult(gid, state, int(row.get('mult_hint') or row.get('mult') or 1))
            r['desired_mult'] = desired
            if r.get('mult') != desired:
                r['geometry_spin_qc'] = 'wrong_spin'
        except RuntimeError:
            r['geometry_spin_qc'] = 'stale_spin_decision'
    r["legacy_g_thermal_eV"] = r.get("g_thermal_eV")
    r["legacy_n_imag"] = r.get("n_imag")
    policy = Path(state_dir) / "thermal_policy.json"
    definition = "legacy_gfn2_gas_rrho"
    selected_protocol = None
    if policy.exists():
        try:
            selected = json.loads(policy.read_text())
            definition = selected["thermal_definition"]
            selected_protocol = selected.get("native_protocol")
        except (OSError, ValueError, KeyError, TypeError):
            definition = "invalid_policy"
    r["thermal_definition"] = definition
    r["thermal_protocol_hash"] = None
    if definition == "legacy_gfn2_gas_rrho":
        r["thermal_qc"] = "imag_modes" if (r.get("n_imag") or 0) > 0 else "ok"
        return r
    r["g_thermal_eV"] = None
    r["n_imag"] = None
    r["thermal_qc"] = "missing_validated_native_thermal"
    if definition != "native_solution_dft_qrrho_spin":
        r["thermal_qc"] = "invalid_thermal_policy"
        return r
    from redox.qm.native_thermal import load_native_thermal, cache_path, NATIVE_THERMAL
    if selected_protocol != NATIVE_THERMAL:
        r["thermal_qc"] = "stale_native_thermal_policy"
        return r
    record = load_native_thermal(state_dir, int(r["charge"]), int(r["mult"]))
    if record is not None:
        r.update(g_thermal_eV=record["g_thermal_total_eV"], n_imag=0,
                 freq_min_cm=min(record["frequencies_cm"]), thermal_qc="ok",
                 thermal_method="solution DFT projected internal-mode qRRHO plus explicit spin partition",
                 thermal_protocol_hash=record["thermal_partition_address"],
                 thermal_partition_protocol=record["thermal_partition_protocol"],
                 thermal_postprocessing_path=record["thermal_postprocessing_path"],
                 native_thermal_record=record)
    return r


def free_energy(r: dict | None):
    """(G_eV, thermal_qc) with G = e_smd + g_thermal from a read_result() dict.

    G is None — the state is INCOMPLETE — if the active-protocol energy or the thermal
    correction is missing. Never substitutes 0 for a missing thermal term (that silently
    mixes electronic-only and free energies in one table). thermal_qc is "ok", "missing",
    or "imag_modes" (xTB RRHO evaluated at a DFT-SMD geometry that is not an xTB stationary
    point had imaginary modes; value kept but flagged for sensitivity review)."""
    if r is None:
        return None, "missing"
    if r.get('geometry_spin_qc', 'ok') != 'ok':
        return None, r['geometry_spin_qc']
    if r.get("thermal_definition") == "native_solution_dft_qrrho_spin" and r.get("thermal_qc") != "ok":
        return None, r.get("thermal_qc", "missing_validated_native_thermal")
    if r.get("thermal_qc") == "invalid_thermal_policy":
        return None, "invalid_thermal_policy"
    e, g = r.get("e_smd_eV"), r.get("g_thermal_eV")
    if g is None:
        return None, "missing"
    qc = "imag_modes" if (r.get("n_imag") or 0) > 0 else "ok"
    if not isinstance(g, (int, float)) or not math.isfinite(g):
        return None, "missing"
    if not isinstance(e, (int, float)) or not math.isfinite(e):
        return None, qc
    return e + g, qc


def group_ids(root: Path = DFT) -> list[str]:
    """Every molecule id with at least one DFT state on disk. Backup dirs (`<id>.bak*`) are
    not molecules and are excluded — the one place this rule lives."""
    return sorted({p.parent.parent.name for p in root.glob("*/*/result.json")
                   if ".bak" not in p.parent.parent.name})


def state_names(gid: str, root: Path = DFT) -> list[str]:
    """Real redox states of `gid`: subdirs holding result.json + opt.xyz. Helper dirs
    (reorg/, sp/, _spincheck*, *.bak*) are never states."""
    d = root / gid
    if not d.exists():
        return []
    return sorted(sd.name for sd in d.iterdir()
                  if sd.is_dir() and not sd.name.startswith(("_", ".")) and ".bak" not in sd.name
                  and sd.name not in ("reorg", "sp")
                  and (sd / "result.json").exists() and (sd / "opt.xyz").exists())


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
        value = float(x)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None

# Compatibility exports for established analysis scripts. New geometry checks belong
# in redox.build.geometry; keeping these names avoids changing historical callers.
from redox.build.geometry import (  # noqa: E402,F401
    _skeleton, smiles_to_xyz_map, geometry_stereo_status, redox_core_xyz_indices,
)
