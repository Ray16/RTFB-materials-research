#!/usr/bin/env python3
"""Compute one checked contact geometry's own electronic and native components."""
import argparse
import gc
import importlib.util
import json
import math
import os
from pathlib import Path
import time


def load(path, name):
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    case, sources = a.plan.parent, a.plan.parent.parent/'sources'
    plan = json.loads(a.plan.read_text())
    digest = load(sources/'hashing.py', 'counterion_thermochemistry_hashes').sha256_file
    def verify():
        if any(digest(path) != h for path, h in plan['source_sha256'].items()):
            raise ValueError('Frozen counterion thermochemistry input changed')
    if (plan['role'] != 'primary_counterion_thermochemistry' or (plan['charge'], plan['mult']) != (1, 1)
            or plan['geometry_repair_permitted'] or plan['native_G_qualified']):
        raise ValueError('Exact unrepaired primary chloride microstate required')
    verify()
    owner = dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in
            ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two disjoint CPUs and matching native pool limits required')
    record = dict(status='counterion_thermochemistry_running', owner=owner, started=time.time(), native_G_qualified=False)
    atomic_write_json(case/'status.json', record)
    try:
        from redox.qm import dft
        from redox.qm.native_thermal import ensure_native_thermal, cache_path, NATIVE_THERMAL
        from redox.core.protocol import ACTIVE_SP
        from redox.validation.protonation import hydrogen_topology_status
        if plan['native_protocol'] != NATIVE_THERMAL or plan['electronic_protocol'] != ACTIVE_SP:
            raise ValueError('Pinned primary electronic/native definitions changed')
        geometry = case/'state/opt.xyz'
        if digest(geometry) != plan['geometry_sha256'] or hydrogen_topology_status(plan['smiles'], geometry) != 'ok':
            raise ValueError('Actual microstate graph or geometry changed')
        helper = load(sources/'helper.py', 'counterion_thermo_screening')
        checkpoints = load(sources/'checkpoints.py', 'counterion_thermo_checkpoints')
        monitor = load(ROOT/'scripts/pipeline/run_pi_stacking.py', 'counterion_thermo_GPU_monitor')
        builder = dft._build_smd_mf
        with helper.corrected_screening(sources/'kernel.cu', case/'screening_implementation'), monitor.MemoryMonitor() as memory:
            dft._build_smd_mf = checkpoints.recording_builder(builder, case/'electronic_SCF_evaluations', owner, save_density=False)
            protocol = plan['electronic_protocol']
            energy = dft.dft_smd(geometry, 1, 1, backend='gpu', do_opt=False, do_freq=False, do_gas=False,
                do_smd=True, sp_xc=protocol['xc'], sp_basis=protocol['basis'], sp_basis_anion=protocol['basis'],
                sp_disp=protocol['disp'], sp_nlc=protocol['nlc'], solvent=protocol['solvent'])
            atomic_write_json(case/'electronic.json', energy)
            if not energy['converged_smd'] or energy['energy_protocol'] != protocol or not math.isfinite(energy['e_smd_eV']):
                raise ValueError('Actual matched electronic energy failed')
            gc.collect(); dft._free_gpu_pool()
            dft._build_smd_mf = checkpoints.recording_builder(builder, case/'SCF_evaluations', owner, save_density=False)
            kind, thermal = ensure_native_thermal(plan['native_id'], plan['native_state'],
                state_dir=case/'state', memory_cap_gib=6)
            if kind != 'computed':
                raise ValueError('Fresh own native thermochemistry required')
            native_path, _ = cache_path(case/'state', 1, 1)
            values = {cutoff: energy['e_smd_eV'] + thermal['nuclear_by_cutoff_eV'][cutoff]
                      + thermal['g_spin_partition_eV'] for cutoff in ('25','50','100')}
            if any(not math.isfinite(value) for value in values.values()):
                raise ValueError('Finite microstate components required')
            atomic_write_json(case/'result.json', dict(key=plan['source_id'], charge=1, mult=1,
                electronic_eV=energy['e_smd_eV'], native_record_path=str(native_path),
                nuclear_by_cutoff_eV=thermal['nuclear_by_cutoff_eV'], spin_partition_eV=thermal['g_spin_partition_eV'],
                G_by_cutoff_eV=values, native_G_qualified=False, association_free_energy_computed=False,
                experimental_populations_inferred=False, material_admission=False))
        dft._build_smd_mf = builder
        if not 0 < memory.peak <= plan['maximum_peak_GPU_MiB']:
            raise ValueError('Prospective electronic/native GPU resource bound exceeded')
        verify()
        record.update(status='counterion_thermochemistry_finished_independent_review_pending', observed_peak_GPU_MiB=memory.peak)
    except BaseException as error:
        record.update(status='counterion_thermochemistry_failed_preserved', error=repr(error))
        raise
    finally:
        record['finished'] = time.time()
        atomic_write_json(case/'status.json', record)


if __name__ == '__main__':
    main()
