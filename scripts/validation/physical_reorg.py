#!/usr/bin/env python
"""Bounded, resumable nonequilibrium PCM validation; never modifies production.

Reference mode compares two independent SCF assemblies on formaldehyde.
Target mode evaluates vertical forward/reverse charge changes at each initial
geometry; no slow-charge transport between different cavities is performed.
Gneq(final at initial geometry; initial solvent) - Geq(final at final geometry)
is a total solution relaxation penalty. Its geometric and solvent pieces are
reported separately; existing gas/outer values must not be added to this total.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import signal
import time

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import DFT, ROOT, HARTREE_EV
from redox.core.protocol import ACTIVE_SP, geom_sha1
from redox.qm.dft import _ecp_for, _kernel_robust
from redox.qm.nonequilibrium_pcm import (
    attach_frozen_slow, independent_hcore_reference, pcm_object, host_array,
    surface_potential,
)

OUT = ROOT / "results/validation/physical_reorg_2026-10-10"
TARGETS = [("anthraquinone_parent", "neu", "red1", "bound_control"),
           ("anthraquinone_parent", "red1", "red2", "gas_unbound_dianion"),
           ("anthraquinone", "neu", "red1", "grafted_conformer_jump")]
REUSE_WORKER_HASH = None


def dense_pcm_operators_on_actual_surface(cm):
    """Independent CPU formulas on the *same* discretized cavity as this event."""
    from pyscf.solvent.pcm import get_D_S, get_F_A
    surface = {k: host_array(v) for k, v in cm.surface.items()}
    D, S = get_D_S(surface, with_S=True, with_D=True)
    _, A = get_F_A(surface)
    f = (float(cm.eps) - 1) / (float(cm.eps) + 1)
    DA = D * A
    K = S - f / (2 * np.pi) * (DA @ S)
    R = -f * (np.eye(len(S)) - DA / (2 * np.pi))
    return K, R


def validate_preserved_solver(path):
    """Approve reuse only when the entire molecular SCF producer is unchanged."""
    old = ast.parse(Path(path).read_text())
    new = ast.parse(Path(__file__).read_text())
    for name in ("make_mf", "solve", "state_info"):
        before = next(n for n in old.body if isinstance(n, ast.FunctionDef) and n.name == name)
        after = next(n for n in new.body if isinstance(n, ast.FunctionDef) and n.name == name)
        if ast.dump(before, include_attributes=False) != ast.dump(after, include_attributes=False):
            raise RuntimeError(f"cannot reuse old SCFs: {name} producer changed")
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, rec):
    temp = path.with_suffix(f".{os.getpid()}.tmp")
    temp.write_text(json.dumps(rec, indent=2, allow_nan=False))
    temp.replace(path)


def make_mf(mol, backend, xc, nlc=None):
    if backend == "gpu":
        from gpu4pyscf import dft
    else:
        from pyscf import dft
    mf = (dft.UKS if mol.spin else dft.RKS)(mol).density_fit()
    mf.xc = xc
    if nlc:
        mf.nlc = nlc
    mf.max_cycle = 200
    mf.conv_tol = 1e-9
    return mf


def solve(mf, dm0=None):
    # Every open-shell case is tested with all three established guesses.
    guesses = ["minao", "atom", "huckel"] if mf.mol.spin else ["minao"]
    records = []
    best = None
    for guess in guesses:
        mf.init_guess = guess
        # gpu4pyscf SCF.kernel silently uses the previous MO density when these
        # are set, defeating an independent init_guess comparison.
        mf.mo_coeff = None
        mf.mo_occ = None
        mf.mo_energy = None
        mf.converged = False
        e, ok = _kernel_robust(mf)
        if ok and (not mf.converged or abs(float(mf.e_tot)-e)>1e-8):
            raise RuntimeError('SCF fallback returned energy without matching density/object')
        rec = dict(guess=guess, e_Ha=e, converged=ok)
        if ok:
            dm = mf.make_rdm1().copy()
            if mf.mol.spin:
                rec["s2"], rec["effective_multiplicity"] = map(float, mf.spin_square())
            if best is None or e < best[0]:
                best = (e, dm, np.asarray(host_array(mf.mo_energy)).tolist(), rec)
        records.append(rec)
    if best is None:
        raise RuntimeError("No converged SCF guess")
    return best, records


def reference(backend="cpu"):
    # Small molecule exercises the full self-consistent field, not only contraction.
    atom = "C 0 0 0; O 0 0 1.21; H .94 0 -.58; H -.94 0 -.58"
    neutral = gto.M(atom=atom, basis="6-31+g*", charge=0, spin=0, verbose=0)
    initial = make_mf(neutral, "cpu", "pbe").PCM(
        pcm_object(neutral, 37.5, "cpu", order=11))
    e0 = float(initial.kernel())
    if not initial.converged:
        raise RuntimeError("Reference initial SCF failed")
    refdm = initial.make_rdm1()
    cation = gto.M(atom=atom, basis="6-31+g*", charge=1, spin=1, verbose=0)
    a = attach_frozen_slow(make_mf(cation, "cpu", "pbe"), refdm,
                           37.5, 1.806, "cpu", order=11)
    ea = float(a.kernel())
    b, constant = independent_hcore_reference(
        make_mf(cation, "cpu", "pbe"), refdm, 37.5, 1.806, order=11)
    eb = float(b.kernel()) + constant
    if not a.converged or not b.converged:
        raise RuntimeError("Independent nonequilibrium reference SCF failed")
    rec = dict(status="passed" if abs(ea-eb)<1e-7 else "failed",
               method="PBE/6-31+G* IEF-PCM", initial_Ha=e0,
               noneq_wrapper_Ha=ea, independent_hcore_Ha=eb,
               difference_eV=(ea-eb)*HARTREE_EV,
               vertical_ionization_eV=(ea-e0)*HARTREE_EV,
               density_max_difference=float(np.max(abs(a.make_rdm1()-b.make_rdm1()))),
               backend="cpu", note="independent energy/SCF assembly, same PCM physics")
    if backend == "gpu":
        c = attach_frozen_slow(make_mf(cation, "gpu", "pbe"), refdm,
                               37.5, 1.806, "gpu", order=11)
        ec = float(c.kernel())
        if not c.converged:
            raise RuntimeError("GPU reference SCF failed")
        rec.update(gpu_noneq_Ha=ec, cpu_gpu_difference_eV=(ec-ea)*HARTREE_EV)
        rec["status"] = "passed" if abs(ec-ea)<1e-6 and rec["status"]=="passed" else "failed"
    atomic_json(OUT / "independent_reference.json", rec)
    if rec["status"] != "passed":
        raise RuntimeError(f"Reference comparison failed: {rec}")
    print(json.dumps(rec), flush=True)


def state_info(gid, state):
    r = json.loads((DFT/gid/state/"result.json").read_text())
    return int(r["charge"]), int(r["mult"])


def cached_scf(gid, state, at_state, backend, reference_state=None, reference_dm=None,
               order=29):
    geom = DFT/gid/at_state/"opt.xyz"
    q, mult = state_info(gid, state)
    key = dict(id=gid, state=state, at_state=at_state, charge=q, mult=mult,
               geom_sha1=geom_sha1(geom), reference_state=reference_state,
               protocol=ACTIVE_SP, model="electrostatic_IEF_PCM", order=order,
               eps_s=37.5, eps_op=1.806, backend=backend,
               code_sha256=hashlib.sha256((ROOT/"src/redox/qm/nonequilibrium_pcm.py").read_bytes()).hexdigest())
    key['worker_code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if reference_dm is not None:
        key["reference_dm_sha256"] = hashlib.sha256(host_array(reference_dm).tobytes()).hexdigest()
    digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:20]
    p = OUT / "scf" / f"{gid}_{state}_at_{at_state}_{reference_state or 'eq'}_{digest}.json"
    p.parent.mkdir(exist_ok=True)
    if p.exists() and p.with_suffix('.npz').exists():
        rec = json.loads(p.read_text())
        if rec.get("key") == key and rec.get("converged") is True:
            dm = np.load(p.with_suffix('.npz'))['dm']
            return rec, dm
    if REUSE_WORKER_HASH:
        # The failed v2 jobs completed their SCFs; only CPU/GPU grid regeneration in
        # the independent contraction failed. Retain those energies and exact densities.
        for previous in sorted(p.parent.glob(f"{gid}_{state}_at_{at_state}_{reference_state or 'eq'}_*.json")):
            old = json.loads(previous.read_text())
            old_key = dict(old.get("key", {}))
            old_hash = old_key.pop("worker_code_sha256", None)
            compare_key = dict(key); compare_key.pop("worker_code_sha256")
            if (old_hash == REUSE_WORKER_HASH and old_key == compare_key
                    and old.get("converged") is True and np.isfinite(old.get("e_Ha", np.nan))
                    and previous.with_suffix('.npz').exists()):
                dm = np.load(previous.with_suffix('.npz'))['dm']
                if not np.isfinite(dm).all():
                    raise RuntimeError(f"nonfinite preserved density {previous}")
                print(f"SCF_REUSE_VERIFIED {previous.name}", flush=True)
                return old, dm
    atoms = read(str(geom))
    symbols = atoms.get_chemical_symbols()
    mol = gto.M(atom=list(zip(symbols, atoms.positions.tolist())),
                basis=ACTIVE_SP["basis"], ecp=_ecp_for(symbols, ACTIVE_SP['basis']),
                charge=q, spin=mult-1, verbose=3, max_memory=32000)
    mf = make_mf(mol, backend, ACTIVE_SP['xc'], ACTIVE_SP['nlc'])
    if reference_dm is None:
        mf = mf.PCM(pcm_object(mol, 37.5, backend, order=order))
    else:
        mf = attach_frozen_slow(mf, reference_dm, 37.5, 1.806, backend, order=order)
    print(f"SCF_START {gid}/{state}@{at_state} ref={reference_state} nao={mol.nao}",flush=True)
    started = time.time()
    best, guesses = solve(mf)
    e, dm, mo, spin = best
    rec = dict(key=key, e_Ha=e, e_eV=e*HARTREE_EV, converged=True,
               guesses=guesses, mo_energy_Ha=mo, s2=spin.get('s2'),
               elapsed_seconds=time.time()-started)
    np.savez_compressed(p.with_suffix('.npz'), dm=host_array(dm))
    atomic_json(p,rec)
    print(f"SCF_COMPLETE {p.name} seconds={rec['elapsed_seconds']:.1f}",flush=True)
    return rec, host_array(dm)


def target(gid, ox, red, role, backend, order):
    entries = {}
    for st in (ox, red):
        entries[st], _ = cached_scf(gid, st, st, backend, order=order)
    rows = []
    for initial, final in ((ox,red),(red,ox)):
        ei, dm0 = cached_scf(gid, initial, initial, backend, order=order)
        efcross, _ = cached_scf(gid, final, initial, backend, order=order)
        eneq, dmneq = cached_scf(gid, final, initial, backend,
                                initial, dm0, order)
        efinal = entries[final]['e_Ha']
        # Native equation must agree with independently contracted density PCM
        # identity at the ACTUAL final nonequilibrium SCF density (not eq density).
        atom = read(str(DFT/gid/initial/'opt.xyz'))
        q,m = state_info(gid,final)
        mol=gto.M(atom=list(zip(atom.get_chemical_symbols(),atom.positions.tolist())),
                  basis=ACTIVE_SP['basis'],charge=q,spin=m-1,verbose=0)
        cm_s=pcm_object(mol,37.5,backend,order=order)
        cm_op=pcm_object(mol,1.806,backend,order=order)
        v0=surface_potential(cm_s,dm0); v1=surface_potential(cm_s,dmneq)
        ne=attach_frozen_slow(make_mf(mol,backend,ACTIVE_SP['xc'],ACTIVE_SP['nlc']),
                              dm0,37.5,1.806,backend,order=order)
        penalty=float(ne.with_solvent._get_vind(dmneq)[0]-cm_s._get_vind(dmneq)[0])
        from redox.properties.lambda_outer import _q_sym
        # Independent CPU dense operators against GPU native iterative solves;
        # surface integrals stay on GPU (CPU 3-center AO blocks are costly).
        # CPU/GPU cavity pruning can differ at nearly buried Lebedev points by one
        # point. Reconstructing a CPU cavity here therefore changes the vector space.
        # Independent NumPy contraction must use the actual event's operators/grid.
        if not np.array_equal(host_array(cm_s.surface['grid_coords']),
                              host_array(cm_op.surface['grid_coords'])):
            raise RuntimeError('static/optical event grids differ')
        Ks, Rs = dense_pcm_operators_on_actual_surface(cm_s)
        Kop, Rop = dense_pcm_operators_on_actual_surface(cm_op)
        dv=host_array(v1-v0)
        if Ks.shape[0] != len(dv) or Kop.shape[0] != len(dv):
            raise RuntimeError('operator and surface potential dimensions differ')
        lam=.5*dv@(_q_sym(Kop,Rop,dv)-_q_sym(Ks,Rs,dv))
        rows.append(dict(initial=initial,final=final,
            total_solution_relaxation_eV=(eneq['e_Ha']-efinal)*HARTREE_EV,
            equilibrium_geometry_component_eV=(efcross['e_Ha']-efinal)*HARTREE_EV,
            selfconsistent_frozen_slow_increment_eV=(eneq['e_Ha']-efcross['e_Ha'])*HARTREE_EV,
            vertical_charge_change_eV=(eneq['e_Ha']-ei['e_Ha'])*HARTREE_EV,
            fixed_density_pekar_penalty_eV=penalty*HARTREE_EV,
            independent_contraction_residual_eV=(penalty-lam)*HARTREE_EV))
    rec=dict(id=gid,couple=f"{ox}->{red}",role=role,method='wb97m-v/def2-tzvpd electrostatic IEF-PCM',
             protocol_distinct_from_smd=True,geometry='existing production SMD geometries',
             validation_status='diagnostic_only',directional=rows,
             self_exchange_total_eV=sum(r['total_solution_relaxation_eV'] for r in rows),
             single_body_symmetric_average_eV=.5*sum(r['total_solution_relaxation_eV'] for r in rows),
             caveat='no independent molecule-matched experimental lambda; no promotion')
    rec['audit_worker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    rec['reused_solver_worker_sha256']=REUSE_WORKER_HASH
    if any(r['total_solution_relaxation_eV']<0 for r in rows):
        rec['physical_flag']='negative_relaxation_nonstationary_endpoint_or_method_problem'
    atomic_json(OUT/f'{gid}_{ox}_{red}.json',rec)
    return rec


def main():
    global REUSE_WORKER_HASH
    p=argparse.ArgumentParser()
    p.add_argument('--backend',choices=['cpu','gpu'],default='cpu')
    p.add_argument('--reference-only',action='store_true')
    p.add_argument('--order',type=int,default=29)
    p.add_argument('--reuse-preserved-solver',type=Path,
                   help='reuse exact method/density records only if SCF producer code matches preserved source')
    a=p.parse_args()
    if a.reuse_preserved_solver:
        REUSE_WORKER_HASH=validate_preserved_solver(a.reuse_preserved_solver)
    OUT.mkdir(parents=True,exist_ok=True)
    def terminate(signum, frame):
        raise KeyboardInterrupt(f'termination signal {signum}')
    signal.signal(signal.SIGTERM,terminate)
    owner=dict(host=os.uname().nodename,pid=os.getpid(),
               process_start=Path('/proc/self/stat').read_text().split()[21],
               cpus=sorted(os.sched_getaffinity(0)),created=time.time())
    atomic_json(OUT/'worker_owner.json',owner)
    reference(a.backend)
    if a.reference_only:
        return
    failures=0
    for gid,ox,red,role in TARGETS:
        claim=OUT/f'{gid}_{ox}_{red}.claim'
        try:
            claim.mkdir()
        except FileExistsError:
            raise RuntimeError(f'Claim exists; inspect owner before restart: {claim}')
        atomic_json(claim/'owner.json',owner)
        try:
            result=target(gid,ox,red,role,a.backend,a.order)
            print(json.dumps(result),flush=True)
            failures=0
        except Exception as exc:
            failures+=1
            atomic_json(OUT/f'{gid}_{ox}_{red}_failure.json',dict(error=repr(exc),owner=owner))
            print(f'TARGET_FAILED {gid}/{ox}->{red}: {exc}',flush=True)
            if failures>=3:
                raise RuntimeError('Three consecutive failures; circuit breaker') from exc
        finally:
            (claim/'owner.json').unlink()
            claim.rmdir()
    atomic_json(OUT/'worker_finished.json',dict(owner=owner,finished=time.time(),
        results_present=[(OUT/f'{g}_{o}_{r}.json').exists() for g,o,r,_ in TARGETS]))
    if not all((OUT/f'{g}_{o}_{r}.json').exists() for g,o,r,_ in TARGETS):
        raise SystemExit(1)


if __name__=='__main__':
    main()
