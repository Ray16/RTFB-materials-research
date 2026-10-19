"""Compute (or adopt) the ACTIVE-protocol single-point energy record for every state.

For each state with an optimized geometry (calcs/dft/<id>/<state>/opt.xyz):
  * record already present          -> skip;
  * result.json was scored at exactly the active protocol on this geometry and RECORDS it
    (xc/basis/nlc/disp/solvent, density fitting, conv_tol, open-shell guess sweep;
    converged)  -> ADOPT its energies (no recompute);
  * otherwise                       -> gas + SMD single points at the active protocol.

Writes calcs/dft/<id>/<state>/sp/<hash>.json (see redox.core.protocol). Never modifies
result.json. Resumable and shardable; run one shard per GPU through gpu_reserve:

  python -m redox.qm.sp --audit                       # what is missing / adoptable
  python -m redox.qm.sp --all --backend gpu --shard 4:0
  python -m redox.qm.sp --only ferrocene --backend gpu
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from redox.core.common import DFT, read_manifest
from redox.core.io import atomic_write_json
from redox.core.protocol import (ACTIVE_SP, geom_sha1, load_record, record_path,
                                 record_provenance_ok, software_versions)

_ENERGY_KEYS = ("e_smd_Ha", "e_smd_eV", "converged_smd", "e_gas_Ha", "e_gas_eV",
                "converged_gas", "gas_homo_eV", "gas_lumo_eV", "gas_s_squared",
                "gas_spin_contam", "anion_unbound", "dG_solv_eV", "density_fit", "conv_tol",
                "scf_guesses", "scf_guess_energies_smd", "scf_guess_energies_gas", "gas_occupation_guesses",
                "smd_occupation_guesses", "smd_homo_eV", "smd_lumo_eV", "smd_s_squared", "smd_spin_contam")


def states(root: Path = DFT, ids=None):
    """[(gid, state, state_dir, charge, mult)] for real redox states (manifest rows with a
    geometry). Ignores helper dirs (reorg/, sp/, _spincheck*, *.bak*)."""
    out = []
    for r in read_manifest():
        gid, st = r["id"], r["state"]
        if ids and gid not in ids and f"{gid}:{st}" not in ids:
            continue
        sd = root / gid / st
        rj = sd / "result.json"
        if not (rj.exists() and (sd / "opt.xyz").exists()):
            continue
        res = json.loads(rj.read_text())
        out.append((gid, st, sd, int(res["charge"]), int(res["mult"])))
    return out


def states_orop():
    """OROP benchmark states (calcs/orop/<sys>/<ox|red>) — same record layout, so the external
    benchmark is scored at exactly the production protocol."""
    root = DFT.parent / "orop"
    out = []
    for sd in sorted(root.glob("*/*/result.json")):
        sd = sd.parent
        if not (sd / "opt.xyz").exists():
            continue
        r = json.loads((sd / "result.json").read_text())
        out.append((f"orop{sd.parent.name}", sd.name, sd, int(r["charge"]), int(r["mult"])))
    return out


def _adoptable(res: dict, proto: dict = ACTIVE_SP, *, charge=None, mult=None,
               geom_hash=None) -> bool:
    """A result.json energy may be adopted only if it PROVES it was computed at the active
    protocol — including density fitting, SCF tolerance and the open-shell guess sweep, which
    result.json files written before 2026-10-02 do not record (anion states optimized before
    RI-J became the default on 2026-09-08 were adopted without it: 7-12 meV off)."""
    if (not isinstance(res, dict) or not geom_hash
            or res.get("energy_geom_sha1") != geom_hash
            or res.get("energy_protocol") != proto
            or type(res.get("charge")) is not int or type(res.get("mult")) is not int
            or res.get("charge") != charge or res.get("mult") != mult
            or type(mult) is not int or mult < 1):
        return False
    if any(type(res.get(key)) not in (int, float) or not math.isfinite(res[key])
           for key in ("e_smd_eV", "e_gas_eV")):
        return False
    guesses = res.get("scf_guesses")
    open_shell = mult > 1
    return (res.get("density_fit") is proto["density_fit"]
            and res.get("conv_tol") == proto["conv_tol"]
            and (not open_shell or (isinstance(guesses, list)
                                   and all(isinstance(guess, str) for guess in guesses)
                                   and len(set(guesses)) > 1))
            and res.get("optimized") is True
            and isinstance(res.get("sp_xc"), str) and res["sp_xc"].lower() == proto["xc"]
            and isinstance(res.get("sp_basis"), str) and res["sp_basis"].lower() == proto["basis"]
            and (res.get("sp_nlc") or None) == proto["nlc"]
            and (res.get("sp_disp") or None) == proto["disp"]
            and res.get("solvent") == proto["solvent"]
            and res.get("converged_smd") is True and res.get("converged_gas") is True
            and res.get("e_smd_eV") is not None and res.get("e_gas_eV") is not None)


def ensure(gid, st, sd: Path, q: int, m: int, backend="gpu", proto: dict = ACTIVE_SP):
    """Make sure the active-protocol record exists. Returns (status, record)."""
    # dft_smd currently fixes these solver settings internally. Never label a
    # calculation with requested settings that the producer cannot apply.
    fixed = ("schema", "density_fit", "conv_tol", "max_cycle", "solvent_model")
    if set(proto) != set(ACTIVE_SP) or any(proto.get(key) != ACTIVE_SP[key] for key in fixed):
        raise ValueError("Unsupported SP solver protocol override")
    rp = record_path(sd, q, m, proto)
    if rp is None:
        raise ValueError(f"{gid}/{st}: missing opt.xyz")
    if rp.exists():
        old = load_record(sd, q, m, proto)
        if old is not None:
            return "present", old
        # unprovable (adopted without density-fitting provenance): keep it, out of the way
        sup = rp.parent / "_superseded"
        sup.mkdir(exist_ok=True)
        rp.rename(sup / f"{rp.stem}.{time.time_ns()}.json")
    source_bytes = (sd / "result.json").read_bytes()
    res = json.loads(source_bytes)
    base = dict(protocol=proto, protocol_hash=rp.stem, geom=str(sd / "opt.xyz"),
                geom_sha1=geom_sha1(sd / "opt.xyz"), charge=q, mult=m, id=gid, state=st)
    if _adoptable(res, proto, charge=q, mult=m, geom_hash=base["geom_sha1"]):
        rec = dict(base, source="adopted:result.json",
                   adopted_result_sha256=hashlib.sha256(source_bytes).hexdigest(),
                   **{k: res.get(k) for k in _ENERGY_KEYS})
        status = "adopted"
    else:
        import redox.qm.dft as D
        t0 = time.time()
        out = D.dft_smd(sd / "opt.xyz", q, m, sp_xc=proto["xc"], sp_basis=proto["basis"],
                        sp_basis_anion=proto["basis"], sp_disp=proto["disp"],
                        sp_nlc=proto["nlc"], solvent=proto["solvent"],
                        do_opt=False, do_gas=True, do_smd=True, do_freq=False,
                        backend=backend)
        rec = dict(base, source=f"computed:{backend}", wall_s=round(time.time() - t0, 1),
                   software=software_versions(),
                   **{k: out.get(k) for k in _ENERGY_KEYS})
        status = "computed"
    rp.parent.mkdir(parents=True, exist_ok=True)
    if geom_sha1(sd / "opt.xyz") != base["geom_sha1"]:
        raise RuntimeError(f"{gid}/{st}: geometry changed during single-point evaluation")
    atomic_write_json(rp, rec)
    return status, rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--only", default=None, help="comma list of gid or gid:state")
    ap.add_argument("--shard", default=None, help="'n:i'")
    ap.add_argument("--backend", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--audit", action="store_true", help="report only; compute nothing")
    ap.add_argument("--orop", action="store_true", help="the OROP benchmark tree instead")
    a = ap.parse_args()
    ids = set(a.only.split(",")) if a.only else None
    todo = states_orop() if a.orop else states(ids=ids)
    if a.orop and ids:
        todo = [t for t in todo if t[0] in ids or f"{t[0]}:{t[1]}" in ids]
    if a.shard:
        n, i = (int(x) for x in a.shard.split(":"))
        todo = [t for k, t in enumerate(todo) if k % n == i]
    if a.audit:
        cnt = {"present": 0, "adoptable": 0, "compute": 0}
        for gid, st, sd, q, m in todo:
            if load_record(sd, q, m) is not None:
                cnt["present"] += 1
            elif _adoptable(json.loads((sd / "result.json").read_text()),
                            charge=q, mult=m, geom_hash=geom_sha1(sd / "opt.xyz")):
                cnt["adoptable"] += 1
            else:
                cnt["compute"] += 1; print(f"  compute {gid}/{st} q={q:+d} m={m}")
        print(f"active protocol {ACTIVE_SP['xc']}/{ACTIVE_SP['basis']}: {cnt}")
        return
    if not (a.all or a.only or a.orop):
        ap.error("use --all, --only, or --audit")
    failed = False
    for gid, st, sd, q, m in todo:
        try:
            status, rec = ensure(gid, st, sd, q, m, backend=a.backend)
            if load_record(sd, q, m) is None:
                failed = True
            print(f"[{status:8s}] {gid}/{st} q={q:+d} m={m} "
                  f"E_smd={rec.get('e_smd_eV')} conv={rec.get('converged_smd')}/"
                  f"{rec.get('converged_gas')}", flush=True)
        except Exception as exc:
            failed = True
            print(f"[FAIL    ] {gid}/{st}: {type(exc).__name__}: {str(exc)[:200]}", flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
