"""Protocol-addressed energy records — the data contract between compute and analysis.

WHY: a state's `result.json` is written once by the geometry optimization and then skipped
forever ("exists -> done"), whatever method, basis, solvent or geometry produced it. Two
failures followed from that:
  * redox ladders mixed basis sets — neutral/cation states at def2-TZVP, anions at
    def2-TZVPD — so the extra diffuse freedom sat on ONE side of every reduction and every
    disproportionation, with no basis-set cancellation;
  * a changed protocol could never be detected, so campaigns could silently mix.

The fix separates GEOMETRY (result.json + opt.xyz, from the optimizer) from ENERGY (single
points at a declared protocol). Every energy record lives at

    calcs/dft/<id>/<state>/sp/<protocol_hash>.json

where the hash covers exactly what determines the number: the coordinates (content hash of
opt.xyz), charge, multiplicity, functional, basis, dispersion/NLC, solvent + solvent model,
density fitting, and SCF convergence settings. Records are never overwritten; a different
protocol or geometry simply addresses a different file. Readers ask for the ACTIVE protocol
and get None if its record is missing or unconverged — they must report INCOMPLETE, never
fall back to an energy computed some other way.

Software versions (pyscf / gpu4pyscf) are stored in each record as metadata but are NOT in
the hash: a version bump should be audited (`python -m redox.qm.sp --audit`), not silently
invalidate every energy.

ACTIVE_SP uses ONE diffuse basis (def2-TZVPD) for EVERY charge state, including the Fc/Fc+
reference, so every redox ladder (and every candidate vs Fc) has a uniform basis.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

SP_SCHEMA = 1

# The single energy protocol every E°, dG_disp and lambda is computed at.
ACTIVE_SP = dict(
    schema=SP_SCHEMA,
    xc="wb97m-v",
    basis="def2-tzvpd",        # uniform: same diffuse basis for all charge states + Fc
    disp=None,
    nlc="vv10",
    solvent="acetonitrile",
    solvent_model="SMD",
    density_fit=True,          # RI-J
    conv_tol=1e-9,             # dft._kernel_robust
    max_cycle=200,
)


def geom_sha1(path: Path) -> str | None:
    try:
        return hashlib.sha1(Path(path).read_bytes()).hexdigest()[:16]
    except OSError:
        return None


def protocol_hash(geom_hash: str, charge: int, mult: int, proto: dict = ACTIVE_SP) -> str:
    key = dict(proto, geom_sha1=geom_hash, charge=int(charge), mult=int(mult))
    blob = json.dumps(key, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:16]


def record_path(state_dir: Path, charge: int, mult: int, proto: dict = ACTIVE_SP):
    """Path of the energy record for the state's CURRENT opt.xyz under `proto` (or None if
    the state has no geometry)."""
    g = geom_sha1(Path(state_dir) / "opt.xyz")
    if g is None:
        return None
    return Path(state_dir) / "sp" / f"{protocol_hash(g, charge, mult, proto)}.json"


def record_provenance_ok(rec: dict, proto: dict = ACTIVE_SP) -> bool:
    """Does the record PROVE it was computed at `proto`? Records computed by redox.qm.sp are
    (they call redox.qm.dft.dft_smd, density-fitted since 2026-09-08). Records ADOPTED from a
    result.json must carry the density-fitting / SCF-tolerance provenance themselves: anion
    states optimized before RI-J became the default were adopted without it (7-12 meV off;
    FINDINGS #24)."""
    if not str(rec.get("source", "")).startswith("adopted"):
        return True
    return (rec.get("density_fit") is proto["density_fit"]
            and rec.get("conv_tol") == proto["conv_tol"])


def load_record(state_dir: Path, charge: int, mult: int, proto: dict = ACTIVE_SP,
                require_converged: bool = True) -> dict | None:
    """The energy record for this state under `proto`, or None if missing / unconverged /
    of unprovable provenance (record_provenance_ok)."""
    p = record_path(state_dir, charge, mult, proto)
    if p is None or not p.exists():
        return None
    rec = json.loads(p.read_text())
    if not record_provenance_ok(rec, proto):
        return None
    if require_converged and not (rec.get("converged_smd") and rec.get("converged_gas")):
        return None
    return rec


def software_versions() -> dict:
    out = {}
    for mod in ("pyscf", "gpu4pyscf"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = None
    return out
