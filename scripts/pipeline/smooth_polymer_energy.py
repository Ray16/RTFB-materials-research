#!/usr/bin/env python
"""Fresh hybrid energies for an independently checked smooth-model geometry.

The alternate model retains its own records. No primary energy cache, thermal
partition, native free energy or candidate-admission certificate is created.
"""
import argparse
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

MODEL = ROOT / 'scripts/pipeline/smooth_cds_control.py'
GC = ROOT / 'scripts/pipeline/smd_gc_workspace.py'


def check_energy(row, task, method_id):
    if (row['charge'] != task['charge'] or row['mult'] != task['mult']
            or row['solvent'] != 'acetonitrile' or row['optimized']
            or row['sp_xc'] != 'wb97m-v' or row['sp_basis'] != 'def2-tzvpd'
            or not row['converged_smd'] or not row['converged_gas']
            or row['experimental_method_id'] != method_id
            or row['primary_candidate_admission'] or row['native_G_computed']):
        raise ValueError('Experimental electronic identity, convergence or scope differs')
    for key in ('e_smd_eV', 'e_gas_eV', 'smd_spin_contam', 'gas_spin_contam'):
        if not math.isfinite(row[key]):
            raise ValueError('Nonfinite experimental electronic energy/spin data')
    if max(abs(row[k]) for k in ('smd_spin_contam', 'gas_spin_contam')) > .1:
        raise ValueError('Experimental electronic spin-contamination gate failed')
    if any(row.get(k) is not None for k in ('G_eV', 'g_thermal_eV', 'G_smd_eV')):
        raise ValueError('An electronic-only calculation cannot supply a free energy')


def prepare(geometry, review, out):
    original = api.read_plan(geometry)
    api.verify_sources(original['producer_sha256']); api.verify_sources(original['input_sha256'])
    proof = json.loads(review.read_text())
    if (proof['status'] != 'alternate_smooth_geometry_and_resource_pilot_independently_checked'
            or proof['campaign'] != str(geometry) or proof['primary_candidate_admission']
            or proof['native_Hessian_or_G_qualified']):
        raise ValueError('A separate completed smooth-model geometry review is required')
    api.verify_sources(proof['source_sha256'])
    model = api.load_module('smooth_energy_declared_model', MODEL)
    identity = hashlib.sha256(json.dumps(model.MODEL, sort_keys=True).encode()).hexdigest()
    if original['experimental_method_id'] != identity or proof['portable_method_id'] != identity:
        raise ValueError('Geometry and energy solvent models must match')
    source = Path(proof['state_dir']) / 'opt.xyz'
    out.mkdir(parents=True, exist_ok=False)
    xyz = out / 'input.xyz'; shutil.copyfile(source, xyz)
    inputs = {**original['input_sha256'], **proof['source_sha256']}
    for p in (geometry / 'plan.json', review, source, xyz):
        inputs[str(p)] = api.digest(p)
    producers = dict(original['producer_sha256'])
    for p in (Path(__file__).resolve(), MODEL, GC, ROOT / 'src/redox/qm/dft.py'):
        producers[str(p)] = api.digest(p)
    plan = dict(original, input_sha256=inputs, producer_sha256=producers,
        raw_root=str(ROOT / 'calcs/polymer_models' / out.parent.name / out.name),
        execution_driver=str(Path(__file__).resolve()), geometry_review=str(review),
        source_geometry_campaign=str(geometry), scored_geometry_sha256=api.digest(xyz),
        maximum_host_RSS_GiB=32, maximum_peak_GPU_MiB=31744,
        prospective_resource_pilot=True, geometry_only_execution=False,
        energy_only_execution=True, native_pipeline_enabled=False, native_G_computed=False,
        purpose='Fresh alternate-model hybrid electronic/resource pilot on exact checked coordinates',
        created=time.time())
    atomic_write_json(out / 'plan.json', plan)
    print(json.dumps(dict(status='experimental_electronic_pilot_prepared', native_G_computed=False)))


@contextmanager
def scf_checkpoints(directory):
    """Checkpoint each gas and SMD SCF factory without changing its equations."""
    from redox.qm import dft
    original = dft._lowest_scf
    count = 0
    def lowest(factory, *args, **kwargs):
        def build():
            nonlocal count
            mf = factory(); folder = directory / f'{count:04d}'; count += 1
            folder.mkdir(parents=True, exist_ok=False)
            mf.chkfile = str(folder / 'SCF.h5')
            atomic_write_json(folder / 'settings.json', dict(xc=mf.xc, basis=repr(mf.mol.basis),
                charge=mf.mol.charge, spin=mf.mol.spin, coordinates_Bohr=mf.mol.atom_coords().tolist(),
                solvent=getattr(getattr(mf, 'with_solvent', None), 'solvent', None), qualified=False))
            previous = getattr(mf, 'callback', None)
            def callback(env):
                if callable(previous): previous(env)
                atomic_write_json(folder / 'progress.json', dict(cycle=env.get('cycle'),
                    E_Ha=float(env['e_tot']), updated=time.time(), qualified=False))
            mf.callback = callback
            return mf
        return original(build, *args, **kwargs)
    dft._lowest_scf = lowest
    try: yield
    finally: dft._lowest_scf = original


def run(campaign):
    from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability
    verify_gate_owner()
    caps = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'XTB_THREADS')
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in caps):
        raise ValueError('Two disjoint CPUs and matching thread caps are required')
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    host_capacity(2, minimum_ram_gib=64)
    gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    guard = json.loads((campaign / 'host_guard/status.json').read_text())
    if guard['pid'] != os.getpid() or guard['maximum_RSS_GiB'] != 32 or guard['minimum_free_GiB'] != 32:
        raise ValueError('The exact same-PID host-memory guard is required')
    task = plan['tasks'][0]
    if len(plan['tasks']) != 1 or (task['charge'], task['mult']) != (0, 1):
        raise ValueError('One neutral singlet control is in scope')
    model = api.load_module('smooth_energy_actual_model', MODEL)
    identity = hashlib.sha256(json.dumps(model.MODEL, sort_keys=True).encode()).hexdigest()
    if plan['experimental_method_id'] != identity or plan['experimental_model'] != model.MODEL:
        raise ValueError('The actual alternate model differs from its declared geometry parent')
    xyz = campaign / 'input.xyz'
    if api.digest(xyz) != plan['scored_geometry_sha256']:
        raise ValueError('Scored coordinates changed after geometry review')
    raw = Path(plan['raw_root']); raw.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(xyz, raw / 'geometry.xyz')
    stage = campaign / 'experimental_energy'; stage.mkdir(exist_ok=False)
    owner = dict(host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    record = dict(status='experimental_hybrid_energy_resource_pilot_running', owner=owner,
        started=time.time(), raw_root=str(raw), native_G_computed=False, primary_candidate_admission=False)
    path = stage / 'status.json'; atomic_write_json(path, record)
    def stop(signum, frame): raise KeyboardInterrupt('Owned experimental energy stop; preserve SCF checkpoints')
    previous = signal.signal(signal.SIGTERM, stop)
    try:
        from redox.qm import dft
        memory = api.load_module('smooth_energy_resource_monitor', ROOT / 'scripts/pipeline/run_pi_stacking.py')
        gc = api.load_module('smooth_energy_GC_context', GC)
        with memory.MemoryMonitor() as used, (stage / 'calculation.log').open('x') as log, redirect_stdout(log), redirect_stderr(log):
            with gc.compact_smd_execution(raw / 'execution', budget_GiB=28), \
                    model.smooth_cds(campaign / 'model', ROOT / 'results/validation/smooth_CDS_2026-10-25/pilot'), \
                    scf_checkpoints(raw / 'SCF_checkpoints'):
                energy = dft.dft_smd(xyz, task['charge'], task['mult'], backend='gpu',
                    do_opt=False, do_freq=False, solvent='acetonitrile')
            energy.update(experimental_model=model.MODEL, experimental_method_id=identity,
                primary_candidate_admission=False, native_G_computed=False)
            check_energy(energy, task, identity)
            atomic_write_json(raw / 'experimental_energy.json', energy)
        if not 0 < used.peak <= plan['maximum_peak_GPU_MiB']:
            raise ValueError('Experimental energy pilot exceeded its GPU resource bound')
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
        record.update(status='experimental_hybrid_energy_finished_independent_review_pending',
            E_smd_eV=energy['e_smd_eV'], E_gas_eV=energy['e_gas_eV'], observed_peak_GPU_MiB=used.peak,
            experimental_method_id=identity, output_sha256={str(p):api.digest(p) for p in
                [raw / 'geometry.xyz', raw / 'experimental_energy.json', campaign / 'model/model.json']})
    except BaseException as error:
        record.update(status='experimental_hybrid_energy_failed_preserved', error=repr(error)); raise
    finally:
        signal.signal(signal.SIGTERM, previous)
        record['finished'] = time.time(); atomic_write_json(path, record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--geometry', type=Path)
    parser.add_argument('--review', type=Path)
    args = parser.parse_args()
    if args.action == 'prepare':
        if not args.geometry or not args.review: parser.error('prepare requires geometry and review')
        prepare(args.geometry.resolve(), args.review.resolve(), args.campaign.resolve())
    else: run(args.campaign.resolve())


if __name__ == '__main__': main()
