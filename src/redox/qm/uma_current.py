"""Current UMA entry point: prefer Nitro, with an explicit V100 legacy fallback.

Frozen replay modules and caches remain untouched. The dispatcher recreates a
fixed-state calculator when composition, charge or multiplicity changes, and
includes electronic state in ASE cache invalidation.
"""
from __future__ import annotations

import argparse
import gc
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm.uma_model_identity import digest

POLICY = ROOT/'config/uma_execution.json'


def select_backend(capability, require_nitro=False):
    major, minor = capability
    if major in (8,9,10,11,12):
        return 'nitro'
    if (major,minor) == (7,0) and not require_nitro:
        return 'legacy'
    raise RuntimeError('No qualified architecture profile; Nitro requires capability 8.x–12.x, with a separate V100 fallback')


def state_aware_calculator(factory, actual_mode, checkpoint_root=None):
    """Keep one predictor in memory and never reuse an ASE cache across states."""
    import numpy as np
    from ase.calculators.calculator import Calculator, all_changes
    from redox.qm.uma_nitro import state_key

    class StateAware(Calculator):
        implemented_properties = ['energy','forces']

        def __init__(self):
            super().__init__()
            self.fixed = None
            self.current_key = None
            self.evaluations = 0
            self.transitions = 0
            self.source_sha256 = digest(__file__)

        def key(self, atoms):
            return state_key([int(z) for z in atoms.numbers], atoms.info.get('charge'), atoms.info.get('spin'))

        def reset(self):
            super().reset()
            if getattr(self,'fixed',None) is not None:
                self.fixed.reset()

        def check_state(self, atoms, tol=1e-15):
            changes = super().check_state(atoms,tol)
            if self.key(atoms) != self.current_key:
                changes.append('electronic_state')
            return changes

        def calculate(self, atoms=None, properties=('energy','forces'), system_changes=all_changes):
            key = self.key(atoms)
            if key != self.current_key:
                self.fixed = None
                gc.collect()
                self.fixed = factory(list(key[0]),key[1],key[2])
                self.current_key = key
                self.transitions += 1
            super().calculate(atoms,properties,system_changes)
            forces = np.asarray(self.fixed.get_forces(atoms),dtype=float)
            energy = float(self.fixed.get_potential_energy(atoms))
            if actual_mode == 'nitro':
                from redox.qm.uma_nitro import check_actual_nitro
                check_actual_nitro(self.fixed.predictor)
            if not math.isfinite(energy) or forces.shape != (len(atoms),3) or not np.isfinite(forces).all():
                raise ValueError('Invalid UMA energy/force result')
            self.results = dict(energy=energy,forces=forces.copy())
            self.evaluations += 1
            if checkpoint_root is not None:
                target=Path(checkpoint_root)/f'state_{self.transitions:04d}'
                target.mkdir(parents=True,exist_ok=True)
                atomic_write_json(target/'latest.json',dict(numbers=list(key[0]),charge=key[1],mult=key[2],
                    coordinates_A=atoms.positions.tolist(),energy_eV=energy,forces_eV_A=forces.tolist(),
                    evaluation=self.evaluations,actual_mode=actual_mode,source_sha256=self.source_sha256,
                    updated=time.time(),optimized_or_qualified=False))

    return StateAware()


def calculator_for_profile(model, device, profile, legacy_factory=None, checkpoint_root=None):
    if profile == 'nitro':
        if device != 'cuda' or model != 'uma-s-1p2p1':
            raise ValueError('The current Nitro profile requires its declared CUDA model')
        from redox.qm.uma_nitro import make_calculator
        return state_aware_calculator(lambda numbers,q,m:make_calculator(numbers,q,m,model=model),'nitro',checkpoint_root)
    if profile != 'legacy':
        raise ValueError('Unknown UMA execution profile')
    if legacy_factory is None:
        from redox.qm.uma import make_calculator as legacy_factory
    return state_aware_calculator(lambda numbers,q,m:legacy_factory(model,device),'legacy',checkpoint_root)


def nitro_state_directory(group, state, policy):
    """Prefer a complete current-model warm start; reject partial/stale data."""
    directory = ROOT/policy['nitro_output']/group/state
    record_path = directory/'result.json'
    if not directory.exists():
        return None
    if not record_path.is_file():
        raise RuntimeError(f'{group}/{state}: partial Nitro warm start is preserved, not replaced by a fallback')
    record = json.loads(record_path.read_text())
    if (record.get('converged') is not True or record.get('id') != group or record.get('state') != state
            or record.get('model_identity',{}).get('checkpoint_sha256') != policy['checkpoint_sha256']
            or record.get('model_identity',{}).get('reference_sha256') != policy['reference_sha256']
            or record.get('UMA_execution',{}).get('actual_mode') != 'nitro'
            or record.get('UMA_execution',{}).get('policy_sha256') != digest(POLICY)
            or digest(directory/'relaxed.xyz') != record.get('geometry_sha256')):
        raise ValueError(f'{group}/{state}: current-model warm-start identity/QC is invalid')
    from redox.core.common import read_manifest
    from redox.qm.uma_nitro import state_key
    from ase.io import read
    rows=[r for r in read_manifest() if r['id']==group and r['state']==state]
    if (len(rows)!=1 or record.get('smiles')!=rows[0]['smiles']
            or record.get('charge')!=int(rows[0]['charge'])):
        raise ValueError('The Nitro warm start belongs to another exact manifest state')
    atoms=read(directory/'relaxed.xyz')
    state_key([int(z) for z in atoms.numbers],record['charge'],record['chosen_mult'])
    from redox.qm.uma import geometry_checks
    if geometry_checks(record['smiles'],directory/'relaxed.xyz',record['charge'],record['chosen_mult']):
        raise ValueError('The current-model warm-start geometry violates its declared graph/state')
    return directory


def run_dft(arguments, policy, require_nitro=False):
    """Use current warm starts while preserving the original DFT spin override."""
    from redox.qm import dft
    previous_geometry, previous_mult = dft.geom_for, dft._uma_mult

    def geometry(group,state):
        selected = nitro_state_directory(group,state,policy)
        if selected is None and require_nitro:
            raise RuntimeError('The required Nitro warm start is missing')
        return selected/'relaxed.xyz' if selected else previous_geometry(group,state)

    def multiplicity(group,state,fallback):
        overrides = ROOT/'config/spin_overrides.json'
        if overrides.exists() and json.loads(overrides.read_text()).get(group,{}).get(state):
            return previous_mult(group,state,fallback)
        selected = nitro_state_directory(group,state,policy)
        return int(json.loads((selected/'result.json').read_text())['chosen_mult']) if selected else previous_mult(group,state,fallback)

    dft.geom_for, dft._uma_mult = geometry, multiplicity
    before = sys.argv
    try:
        sys.argv = ['redox.qm.dft',*arguments]
        dft.main()
    finally:
        sys.argv = before
        dft.geom_for, dft._uma_mult = previous_geometry, previous_mult


def run_uma(arguments, profile, policy):
    from redox.qm import uma
    original_factory, original_identity, original_scan, original_out = uma.make_calculator, uma.model_identity, uma.relax_state_scan, uma.OUT
    calculators = []

    def factory(model,device):
        calc = calculator_for_profile(model,device,profile,legacy_factory=original_factory,
            checkpoint_root=checkpoint_parent/f'calculator_{len(calculators):04d}')
        calculators.append(calc)
        return calc

    def identity(model):
        if profile == 'legacy':
            return original_identity(model)
        from redox.qm.uma_model_identity import model_identity
        value = model_identity(model)
        if (value['checkpoint_sha256'] != policy['checkpoint_sha256']
                or value['reference_sha256'] != policy['reference_sha256']):
            raise ValueError('The configured Nitro weights or atomic references differ')
        value['execution_profile'] = dict(actual_mode='nitro',semantics='official_current_UMA_S1.2.1',policy_sha256=digest(POLICY))
        return value

    def scan(calc,*args,**kwargs):
        result = original_scan(calc,*args,**kwargs)
        result['UMA_execution'] = dict(requested_mode='nitro',actual_mode=profile,
            state_transitions=calc.transitions,evaluations=calc.evaluations,
            source_sha256=digest(__file__),policy_sha256=digest(POLICY),
            interpretation='Warm-start model output; not a DFT energy or solution free energy')
        return result

    uma.make_calculator,uma.model_identity,uma.relax_state_scan = factory,identity,scan
    uma.OUT = ROOT/policy['nitro_output'] if profile == 'nitro' else ROOT/policy['legacy_output']
    receipts = ROOT/'calcs/uma_profiles/run_receipts';receipts.mkdir(parents=True,exist_ok=True)
    receipt = receipts/f'{time.time_ns()}_{os.getpid()}.json'
    checkpoint_parent = ROOT/'calcs/uma_profiles/checkpoints'/receipt.stem
    report = dict(started=time.time(),requested_backend='nitro',actual_backend=profile,
        arguments=arguments,output_root=str(uma.OUT),source_sha256=digest(__file__),policy_sha256=digest(POLICY))
    atomic_write_json(receipt,report)
    before = sys.argv
    try:
        sys.argv = ['redox.qm.uma',*arguments]
        uma.main()
        report['status']='UMA_command_returned'
    except BaseException as error:
        report.update(status='UMA_command_failed_preserved',error=repr(error));raise
    finally:
        report.update(finished=time.time(),actual_inference_evaluations=sum(c.evaluations for c in calculators))
        atomic_write_json(receipt,report)
        sys.argv = before
        uma.make_calculator,uma.model_identity,uma.relax_state_scan,uma.OUT = original_factory,original_identity,original_scan,original_out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('uma','dft'),default='uma')
    parser.add_argument('--require-nitro',action='store_true')
    parser.add_argument('--show-policy',action='store_true')
    args,remaining=parser.parse_known_args()
    policy=json.loads(POLICY.read_text())
    if args.show_policy:
        print(json.dumps(policy,indent=2));return
    if args.stage == 'dft':
        python=policy['legacy_python']
        if Path(sys.prefix).resolve() != Path(python).parent.parent.resolve():
            environment=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
            os.execve(python,[python,'-m','redox.qm.uma_current',*sys.argv[1:]],environment)
        backend_index=remaining.index('--backend') if '--backend' in remaining else -1
        if '--backend=gpu' in remaining or (0 <= backend_index < len(remaining)-1 and remaining[backend_index+1]=='gpu'):
            from redox.workflows.resources import verify_gate_owner
            verify_gate_owner()
        run_dft(remaining,policy,args.require_nitro);return
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    gpu=os.environ['CUDA_VISIBLE_DEVICES']
    raw=subprocess.check_output(['nvidia-smi','--id='+gpu,'--query-gpu=compute_cap','--format=csv,noheader,nounits'],text=True,timeout=15).strip()
    major,minor=map(int,raw.split('.'))
    profile=select_backend((major,minor),args.require_nitro)
    python=policy['nitro_python'] if profile == 'nitro' else policy['legacy_python']
    wanted_prefix=Path(python).parent.parent
    needs_library_path = profile == 'nitro' and not os.environ.get('LD_LIBRARY_PATH','').startswith(policy['nitro_library_path'])
    if Path(sys.prefix).resolve() != wanted_prefix.resolve() or needs_library_path:
        environment=dict(os.environ,PYTHONPATH=str(ROOT/'src'),FAIRCHEM_CACHE_DIR=policy['fairchem_cache'])
        if profile == 'nitro':environment['LD_LIBRARY_PATH']=policy['nitro_library_path']+os.pathsep+environment.get('LD_LIBRARY_PATH','')
        os.execve(python,[python,'-m','redox.qm.uma_current',*sys.argv[1:]],environment)
    os.environ['FAIRCHEM_CACHE_DIR']=policy['fairchem_cache']
    print(f'[UMA] preferred=nitro actual={profile}',flush=True)
    run_uma(remaining,profile,policy)


if __name__=='__main__':main()
