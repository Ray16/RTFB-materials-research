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
  python -m redox.qm.uma --model uma-s-1p2p1 --device cuda            # all states
  python -m redox.qm.uma --only pyridinium                           # one group
  python -m redox.qm.uma --only pyridinium:red                       # one state
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import tempfile

from redox.core.common import LIBRARY, read_manifest
from redox.core.common import UMA as OUT
# uma-s-1p2p1 (newest UMA, v1.2.1) works in fairchem 2.21 via REGISTRATION: its checkpoint
# is the same architecture as uma-s-1p2 (patched weights, ~0.5 kJ/mol apart), so we add a
# registry entry pointing to uma-s-1p2p1.pt and load it through uma-s-1p2's compatible
# config. (Direct path-load fails because the checkpoint's embedded config has newer
# HydraModel kwargs.) ensure_registered() does this idempotently. torch 2.8+cu128 runs on
# the 12.4 driver via CUDA minor-version compat.
DEFAULT_MODEL = "uma-s-1p2p1"
CACHE_SCHEMA = 1


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def model_identity(model):
    """Fingerprint locally cached weights, atomic references and loader version.

    No checkpoint download or model/GPU initialization is needed for a resume check.
    The same registry entry and HF cache used by the installed loader are inspected.
    """
    from huggingface_hub import hf_hub_download
    dist = importlib.metadata.distribution('fairchem-core')
    registry = Path(dist.locate_file('fairchem/core/calculate/pretrained_models.json'))
    entries = json.loads(registry.read_text())
    entry = entries.get(model)
    if entry is None and model == DEFAULT_MODEL and 'uma-s-1p2' in entries:
        entry = dict(entries['uma-s-1p2'], filename='uma-s-1p2p1.pt')
    if entry is None:
        raise RuntimeError(f'{model}: model registry entry is unavailable')
    cache = os.environ.get('FAIRCHEM_CACHE_DIR', str(Path.home() / '.cache/fairchem'))
    def local_file(record):
        return Path(hf_hub_download(repo_id=entry['repo_id'], filename=record['filename'],
                    subfolder=record.get('subfolder'), revision=entry.get('revision'),
                    cache_dir=cache, local_files_only=True))
    checkpoint = local_file(entry)
    references = {key: file_sha256(local_file(entry[key]))
                  for key in ('atom_refs', 'form_elem_refs') if entry.get(key)}
    return dict(model=model, checkpoint_sha256=file_sha256(checkpoint),
                reference_sha256=references, loader_version=dist.version,
                registry_entry=entry, task_name='omol')


def seed_identity(conf_paths):
    """Ordered ensemble identity; changing any seed invalidates automatic reuse."""
    return [{'name': Path(p).name, 'sha256': file_sha256(p)} for p in conf_paths]


def geometry_checks(smiles, path, charge, mult):
    """Independent input-graph/composition/stereo checks on an optimized XYZ."""
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from redox.core.common import _skeleton, geometry_stereo_status
    lines = Path(path).read_text().splitlines()
    plain = '\n'.join(lines[:2] + [' '.join(line.split()[:4]) for line in lines[2:]]) + '\n'
    actual = Chem.MolFromXYZBlock(plain)
    declared = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if actual is None or declared is None:
        return {'geometry_unreadable': True}
    rdDetermineBonds.DetermineConnectivity(actual, charge=0)
    expected = _skeleton(Chem.RemoveHs(declared))
    observed = _skeleton(Chem.RemoveHs(actual, sanitize=False))
    electrons = sum(a.GetAtomicNum() for a in actual.GetAtoms()) - charge
    errors = dict(composition_mismatch=sorted(a.GetAtomicNum() for a in actual.GetAtoms()) !=
                      sorted(a.GetAtomicNum() for a in declared.GetAtoms()),
                  connectivity_mismatch=not (observed.HasSubstructMatch(expected) and expected.HasSubstructMatch(observed)),
                  spin_parity_mismatch=mult < 1 or electrons < mult - 1 or (electrons - (mult - 1)) % 2 != 0)
    stereo = geometry_stereo_status(smiles, path)
    errors['stereochemistry_mismatch'] = stereo not in ('ok', 'not_applicable')
    return {key: value for key, value in errors.items() if value}


def cache_errors(record, row, xyz, conf_paths, identity, requested_fmax):
    """Return explicit reasons an existing warm start cannot be reused."""
    hint = int(row.get('mult_hint') or row.get('mult') or 1)
    wanted = dict(cache_schema=CACHE_SCHEMA, model_identity=identity, model=identity['model'],
                  seed_identity=seed_identity(conf_paths), smiles=row['smiles'],
                  charge=int(row['charge']), mult_hint=hint, n_e=int(row['n_e']),
                  id=row['id'], state=row['state'])
    errors = [key + '_mismatch_or_missing' for key, value in wanted.items() if record.get(key) != value]
    for key in ('energy_eV', 'fmax'):
        if type(record.get(key)) not in (int, float) or not math.isfinite(record[key]):
            errors.append(key + '_invalid')
    if record.get('converged') is not True or not isinstance(record.get('fmax'), (int, float)) or record['fmax'] > requested_fmax:
        errors.append('convergence_requirement_not_met')
    if not Path(xyz).exists():
        return errors + ['missing_geometry']
    if record.get('geometry_sha256') != file_sha256(xyz):
        errors.append('geometry_changed_or_unfingerprinted')
    from rdkit import Chem
    mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
    electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - int(row['charge'])
    mults = candidate_mults(electrons, hint)
    if record.get('n_elec') != electrons or record.get('mults_scanned') != mults or record.get('chosen_mult') not in mults:
        errors.append('electron_or_spin_scan_identity_mismatch')
    if record.get('n_conf') != len(conf_paths) or type(record.get('conf_idx')) is not int or not 0 <= record['conf_idx'] < len(conf_paths):
        errors.append('conformer_selection_mismatch')
    if not errors:
        errors.extend(geometry_checks(row['smiles'], xyz, int(row['charge']), record['chosen_mult']))
    return errors


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
    tmp = reg.with_suffix(".json.tmp")   # atomic write (avoids parallel-shard races)
    tmp.write_text(json.dumps(d, indent=4))
    import os
    os.replace(tmp, reg)


def make_calculator(model: str, device: str):
    """Build a FAIRChem OMol calculator. API verified against the installed fairchem-core
    before first run (see archive/scripts/probe_uma.py)."""
    from fairchem.core import pretrained_mlip, FAIRChemCalculator
    ensure_registered(model)
    predictor = pretrained_mlip.get_predict_unit(model, device=device)
    return FAIRChemCalculator(predictor, task_name="omol")


def candidate_mults(n_electrons: int, hint: int | None = None):
    """Spin multiplicities to scan for a given electron count. Even-electron species get
    singlet+triplet (the real ground-state question, e.g. reduced dications/dianions);
    odd-electron get doublet. The config hint is always included."""
    base = [1, 3] if n_electrons % 2 == 0 else [2]
    if hint and hint not in base:
        base.append(hint)
    return sorted(set(base))


def relax_one(calc, atoms, charge: int, mult: int, fmax: float, steps: int):
    """Relax one geometry under (charge, mult). Returns (energy_eV, fmax, n_steps, atoms)."""
    from ase.optimize import BFGS
    a = atoms.copy()
    a.info["charge"] = int(charge)   # OMol conditioning: total charge + spin multiplicity
    a.info["spin"] = int(mult)
    a.calc = calc
    opt = BFGS(a, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    f = a.get_forces()
    fmax_final = float((f ** 2).sum(axis=1).max() ** 0.5)
    return float(a.get_potential_energy()), fmax_final, int(opt.get_number_of_steps()), a


def relax_state_scan(calc, conf_paths, charge: int, mult_hint: int, fmax: float, steps: int):
    """Scan every (conformer x candidate multiplicity), relax each with UMA, and return the
    global lowest-energy result + diagnostics (chosen mult, singlet-triplet gap, conformer
    spread) so the spin ground state and best conformer are DETERMINED, not assumed."""
    from ase.io import read
    confs = [read(str(p)) for p in conf_paths]
    n_elec = int(sum(confs[0].get_atomic_numbers())) - int(charge)
    mults = candidate_mults(n_elec, mult_hint)

    best = None                      # (energy, mult, conf_idx, atoms, fmax, nsteps)
    per_mult_best = {}               # mult -> best energy (for spin gap)
    for m in mults:
        for ci, atoms in enumerate(confs):
            e, fm, ns, a = relax_one(calc, atoms, charge, m, fmax, steps)
            if not math.isfinite(e) or not math.isfinite(fm) or fm > fmax:
                continue  # an unconverged low energy cannot win the spin/conformer scan
            if m not in per_mult_best or e < per_mult_best[m]:
                per_mult_best[m] = e
            if best is None or e < best[0]:
                best = (e, m, ci, a, fm, ns)
    if best is None or set(per_mult_best) != set(mults):
        raise RuntimeError('UMA scan lacks a converged conformer for every requested multiplicity')
    e, m, ci, a, fm, ns = best
    others = sorted(v for k, v in per_mult_best.items() if k != m)
    spin_gap = round(others[0] - e, 4) if others else None
    return dict(energy_eV=e, chosen_mult=m, mult_hint=mult_hint, conf_idx=ci,
                fmax=fm, n_steps=ns, converged=bool(fm <= fmax),
                charge=int(charge), n_elec=n_elec, mults_scanned=mults,
                spin_gap_eV=spin_gap, per_mult_eV={k: round(v, 4) for k, v in per_mult_best.items()},
                n_conf=len(confs), atoms=a)


def write_xyz(atoms, path: Path, comment: str):
    from ase.io import write
    path.parent.mkdir(parents=True, exist_ok=True)
    write(str(path), atoms, format='xyz', comment=comment)


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
    if not math.isfinite(args.fmax) or args.fmax <= 0 or args.steps < 1:
        ap.error('fmax must be finite and positive; steps must be positive')

    rows = read_manifest()
    if args.only:
        gid, _, st = args.only.partition(":")
        rows = [r for r in rows if r["id"] == gid and (not st or r["state"] == st)]
        if not rows:
            raise SystemExit(f"no states match --only {args.only!r}")
    if args.shard:
        n, i = (int(x) for x in args.shard.split(":"))
        if n < 1 or not 0 <= i < n:
            ap.error('shard requires n>=1 and 0<=i<n')
        rows = [r for k, r in enumerate(rows) if k % n == i]

    calc = None  # resume fingerprints local files without initializing the model
    identity = None
    for r in rows:
        gid, st = r["id"], r["state"]
        outdir = OUT / gid / st
        res_json = outdir / "result.json"
        conf_paths = sorted((LIBRARY / gid / "conformers").glob("conf_*.xyz"))
        if not conf_paths:
            raise RuntimeError(f"{gid}/{st}: no conformer ensemble (run the candidate builder)")
        if res_json.exists() and not args.force:
            if identity is None:
                identity = model_identity(args.model)
            old = json.loads(res_json.read_text())
            errors = cache_errors(old, r, outdir / 'relaxed.xyz', conf_paths, identity, args.fmax)
            if errors:
                raise RuntimeError(f"{gid}/{st}: UMA cache is unverifiable/stale ({'; '.join(errors)}); "
                                   'preserved for review, not automatically replaced')
            print(f"[skip] {gid}/{st} (done)")
            continue
        if calc is None:
            print(f"[load] {args.model} on {args.device}")
            calc = make_calculator(args.model, args.device)
            identity = model_identity(args.model)
        mult_hint = int(r.get("mult_hint") or r.get("mult") or 1)
        print(f"[run ] {gid}/{st} q={r['charge']} scan {len(conf_paths)} confs x mults ...",
              flush=True)
        res = relax_state_scan(calc, conf_paths, int(r["charge"]), mult_hint,
                               args.fmax, args.steps)
        atoms = res.pop("atoms")
        write_xyz(atoms, outdir / "relaxed.xyz",
                  comment=f"{gid}/{st} q={r['charge']} m={res['chosen_mult']} "
                          f"E={res['energy_eV']:.6f}eV conf={res['conf_idx']} model={args.model}")
        problems = geometry_checks(r['smiles'], outdir / 'relaxed.xyz', int(r['charge']), res['chosen_mult'])
        if problems:
            raise RuntimeError(f'{gid}/{st}: optimized UMA geometry violates input identity: {problems}')
        res.update(id=gid, state=st, model=args.model,
                   n_e=int(r["n_e"]), smiles=r["smiles"], cache_schema=CACHE_SCHEMA,
                   model_identity=identity, seed_identity=seed_identity(conf_paths),
                   geometry_sha256=file_sha256(outdir / 'relaxed.xyz'),
                   requested_fmax=args.fmax, optimizer='ASE BFGS', optimizer_step_limit=args.steps)
        with tempfile.NamedTemporaryFile(mode='w',dir=outdir,prefix='.result-',suffix='.tmp',delete=False) as stream:
            staged=Path(stream.name)
            json.dump(res,stream,indent=2)
        try:
            staged.replace(res_json)
        finally:
            staged.unlink(missing_ok=True)
        flag = " <-- mult != hint" if res["chosen_mult"] != mult_hint else ""
        print(f"[done] {gid}/{st} E={res['energy_eV']:.4f} mult={res['chosen_mult']} "
              f"(hint {mult_hint}){flag} S-T gap={res['spin_gap_eV']} conf={res['conf_idx']} "
              f"conv={res['converged']}")
        if not res["converged"]:
            raise RuntimeError(f"{gid}/{st}: UMA relaxation did not converge")


if __name__ == "__main__":
    main()
