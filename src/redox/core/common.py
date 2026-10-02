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
               "gas_spin_contam", "anion_unbound", "dG_solv_eV")


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
    e, g = r.get("e_smd_eV"), r.get("g_thermal_eV")
    if g is None:
        return None, "missing"
    if e is None:
        return None, "ok"
    qc = "imag_modes" if (r.get("n_imag") or 0) > 0 else "ok"
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
        return float(x)
    except (TypeError, ValueError):
        return None

# --- SMILES <-> stored-geometry atom mapping -------------------------------------------------
# WHY THIS EXISTS: build_candidates.py embeds the mol returned by decorate() but records
# Chem.MolToSmiles(mol) -- the CANONICAL SMILES -- in the manifest. Canonicalization reorders
# atoms, so for every GRAFTED candidate the manifest-SMILES atom indices do NOT correspond to
# the atom order in library/*.xyz or calcs/dft/*/opt.xyz. Standalone entries (built from a
# direct SMILES) happen to round-trip to identity, so anything that indexes a geometry by
# SMILES index silently works on the standalones and is scrambled on the grafted molecules --
# which are the actual screening candidates. Always go through this mapper.
#
# Bond ORDERS cannot be perceived reliably from geometry for charged aromatics (RDKit's
# DetermineBonds fails on viologen dications and radical anions), but CONNECTIVITY can, so the
# match is done on the element+connectivity skeleton with charges and bond orders stripped.

def _skeleton(mol):
    """Element + connectivity only: no formal charges, no explicit-H counts, no bond orders."""
    from rdkit import Chem
    rw = Chem.RWMol(mol)
    for a in rw.GetAtoms():
        a.SetFormalCharge(0); a.SetNumExplicitHs(0); a.SetNoImplicit(True); a.SetIsAromatic(False)
    for b in rw.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
    return rw.GetMol()


def smiles_to_xyz_map(smiles: str, xyz_path):
    """Heavy-atom index in `smiles` -> heavy-atom index in the geometry at `xyz_path`.

    Returns a tuple m where m[i] is the xyz index of SMILES heavy atom i, or None if the
    geometry cannot be read / matched. Verified to map all 12 current candidate molecules.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    try:
        ref = Chem.MolFromXYZFile(str(xyz_path))
        if ref is None:
            return None
        rdDetermineBonds.DetermineConnectivity(ref, charge=0)
        ref = _skeleton(Chem.RemoveHs(ref, sanitize=False))
        q = Chem.MolFromSmiles(smiles)
        if q is None:
            return None
        return ref.GetSubstructMatch(_skeleton(Chem.RemoveHs(q)), useChirality=False) or None
    except Exception:
        return None


def redox_core_xyz_indices(smiles: str, xyz_path, tether_smarts="[CH3]c1ccc([CH2])cc1"):
    """xyz heavy-atom indices of the REDOX CORE, i.e. ring atoms outside the Merrifield
    benzyl tether plus the exocyclic atoms double-bonded to them (quinone / imide C=O).

    Acyclic substituents (N-methyl, N-ethyl, ammoniopropyl, methoxy) are EXCLUDED: they are
    conformational degrees of freedom, not redox-active, and letting them into an RMSD is
    what makes a whole-molecule test too blunt to separate a tether rotation from a genuine
    inner-sphere distortion. Returns (core_idx, tether_idx) as sorted lists of xyz indices.
    """
    from rdkit import Chem
    m = smiles_to_xyz_map(smiles, xyz_path)
    if m is None:
        return None, None
    q = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
    teth = set(q.GetSubstructMatch(Chem.MolFromSmarts(tether_smarts)))
    ri = q.GetRingInfo()
    ring = {i for i in range(q.GetNumAtoms()) if ri.NumAtomRings(i) > 0}
    core = ring - teth
    for a in list(core):
        for b in q.GetAtomWithIdx(a).GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE:
                o = b.GetOtherAtomIdx(a)
                if o not in ring:
                    core.add(o)
    return sorted(m[i] for i in core), sorted(m[i] for i in teth)
