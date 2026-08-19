"""UMA (FAIRChem / OMol) gas-phase geometry optimization for each redox state.

Reads library/manifest.csv, and for every (id, state) relaxes the seed geometry with UMA
under that state's explicit charge + spin multiplicity, writing:

  calcs/uma/<id>/<state>/relaxed.xyz    optimized geometry
  calcs/uma/<id>/<state>/result.json    energy (eV), fmax, n_steps, converged, charge, mult

Resumable: a state whose result.json exists is skipped. GPU is selected via the --device
flag or CUDA_VISIBLE_DEVICES so states can be fanned across GPUs by an external launcher.

NOTE: UMA/OMol is gas-phase. Charged states are pre-optimized here only as a warm start;
solvated geometries come from the DFT+SMD step (see docs/MODELING.md).

Usage:
  python -m redox.uma --model uma-s-1p2p1 --device cuda            # all states
  python -m redox.uma --only pyridinium                           # one group
  python -m redox.uma --only pyridinium:red                       # one state
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIBRARY = ROOT / "library"
OUT = ROOT / "calcs" / "uma"
# uma-s-1p2p1 (newest UMA, v1.2.1) works in fairchem 2.21 via REGISTRATION: its checkpoint
# is the same architecture as uma-s-1p2 (patched weights, ~0.5 kJ/mol apart), so we add a
# registry entry pointing to uma-s-1p2p1.pt and load it through uma-s-1p2's compatible
# config. (Direct path-load fails because the checkpoint's embedded config has newer
# HydraModel kwargs.) ensure_registered() does this idempotently. torch 2.8+cu128 runs on
# the 12.4 driver via CUDA minor-version compat.
DEFAULT_MODEL = "uma-s-1p2p1"


def ensure_registered(model: str):
    """Idempotently register uma-s-1p2p1 in fairchem's model registry by cloning the
    uma-s-1p2 entry with the 1p2p1 checkpoint filename. No-op for other models."""
    if model != "uma-s-1p2p1":
        return
    import json
    import fairchem.core.calculate as fcc
    reg = Path(fcc.__file__).parent / "pretrained_models.json"
    d = json.loads(reg.read_text())
    if "uma-s-1p2p1" in d:
        return
    if "uma-s-1p2" not in d:
        raise RuntimeError("uma-s-1p2 not in fairchem registry; cannot derive 1p2p1")
    entry = dict(d["uma-s-1p2"]); entry["filename"] = "uma-s-1p2p1.pt"
    d["uma-s-1p2p1"] = entry
    reg.write_text(json.dumps(d, indent=4))


def read_manifest():
    with (LIBRARY / "manifest.csv").open() as f:
        return list(csv.DictReader(f))


def make_calculator(model: str, device: str):
    """Build a FAIRChem OMol calculator. API verified against the installed fairchem-core
    before first run (see scripts/probe_uma.py)."""
    from fairchem.core import pretrained_mlip, FAIRChemCalculator
    ensure_registered(model)
    predictor = pretrained_mlip.get_predict_unit(model, device=device)
    return FAIRChemCalculator(predictor, task_name="omol")


def relax_state(calc, xyz_path: Path, charge: int, mult: int, fmax: float, steps: int):
    """Relax one geometry with UMA under (charge, mult). Returns dict of results."""
    from ase.io import read
    from ase.optimize import BFGS
    atoms = read(str(xyz_path))
    # OMol conditioning: total charge + spin multiplicity (2S+1) carried on atoms.info.
    atoms.info["charge"] = int(charge)
    atoms.info["spin"] = int(mult)
    atoms.calc = calc
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    f = atoms.get_forces()
    fmax_final = float((f ** 2).sum(axis=1).max() ** 0.5)
    return dict(
        energy_eV=float(atoms.get_potential_energy()),
        fmax=fmax_final,
        n_steps=int(opt.get_number_of_steps()),
        converged=bool(fmax_final <= fmax),
        charge=int(charge), mult=int(mult),
        atoms=atoms,
    )


def write_xyz(atoms, path: Path, comment: str):
    from ase.io import write
    path.parent.mkdir(parents=True, exist_ok=True)
    write(str(path), atoms, comment=comment)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--fmax", type=float, default=0.02)
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--only", default=None,
                    help="restrict to 'group' or 'group:state'")
    ap.add_argument("--shard", default=None,
                    help="'n:i' — process only rows where index%%n==i (for GPU fan-out)")
    ap.add_argument("--force", action="store_true", help="recompute even if result exists")
    args = ap.parse_args()

    rows = read_manifest()
    if args.only:
        gid, _, st = args.only.partition(":")
        rows = [r for r in rows if r["id"] == gid and (not st or r["state"] == st)]
        if not rows:
            raise SystemExit(f"no states match --only {args.only!r}")
    if args.shard:
        n, i = (int(x) for x in args.shard.split(":"))
        rows = [r for k, r in enumerate(rows) if k % n == i]

    calc = None  # lazily built so --only/resume listing is cheap and import errors surface late
    for r in rows:
        gid, st = r["id"], r["state"]
        outdir = OUT / gid / st
        res_json = outdir / "result.json"
        if res_json.exists() and not args.force:
            print(f"[skip] {gid}/{st} (done)")
            continue
        seed = LIBRARY / gid / f"{st}.xyz"
        if calc is None:
            print(f"[load] {args.model} on {args.device}")
            calc = make_calculator(args.model, args.device)
        print(f"[run ] {gid}/{st} q={r['charge']} m={r['mult']} ...", flush=True)
        res = relax_state(calc, seed, int(r["charge"]), int(r["mult"]),
                          args.fmax, args.steps)
        atoms = res.pop("atoms")
        write_xyz(atoms, outdir / "relaxed.xyz",
                  comment=f"{gid}/{st} q={r['charge']} m={r['mult']} "
                          f"E={res['energy_eV']:.6f}eV model={args.model}")
        res.update(id=gid, state=st, model=args.model,
                   n_e=int(r["n_e"]), smiles=r["smiles"])
        res_json.write_text(json.dumps(res, indent=2))
        print(f"[done] {gid}/{st} E={res['energy_eV']:.4f} eV "
              f"fmax={res['fmax']:.4f} steps={res['n_steps']} conv={res['converged']}")


if __name__ == "__main__":
    main()
