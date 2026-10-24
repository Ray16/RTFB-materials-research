#!/usr/bin/env python
"""Qualify CPU PCM tile contractions and run a separately bounded large point."""
import argparse
from contextlib import nullcontext
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

ADAPTER = ROOT/'scripts/pipeline/cpu_pcm_memory.py'


def controls(source, out):
    with allocation_guard():
        host_capacity(2, minimum_ram_gib=16)
        os.sched_setaffinity(0,cpu_allocation(2))
    import numpy as np
    import pyscf.lib
    from pyscf import scf
    from redox.qm import dft
    previous = json.loads((source/'status.json').read_text())
    if previous['status'] != 'two_spin_CPU_GPU_energy_and_gradient_equivalence_passed':
        raise ValueError('Both original CPU/GPU controls required')
    api.verify_sources(previous['source_sha256'])
    out.mkdir(parents=True,exist_ok=False)
    adapter = api.load_module('CPU_PCM_tile_contract_controls',ADAPTER)
    report = dict(status='CPU_PCM_tile_controls_running', started=time.time(),
        host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        maximum_matrix_error_tolerance=1e-10,cases=[],source_sha256={str(p):api.digest(p) for p in
            (Path(__file__).resolve(),ADAPTER,source/'status.json')})
    atomic_write_json(out/'status.json',report)
    pyscf.lib.num_threads(2)
    try:
        for mult in (1,2):
            checkpoint = source/f'geometry_mult{mult}.h5'
            mol, saved = scf.chkfile.load_scf(str(checkpoint))
            dm = scf.uhf.make_rdm1(saved['mo_coeff'],saved['mo_occ']).sum(axis=0) if mult == 2 else scf.hf.make_rdm1(saved['mo_coeff'],saved['mo_occ'])
            outputs = []
            for mode in ('original','bounded'):
                context = adapter.bounded_pcm_memory(out/f'mult{mult}_{mode}',.001) if mode == 'bounded' else nullcontext()
                with context:
                    mf = dft._build_smd_mf(mol,'r2scan','acetonitrile','cpu',disp='d4')
                    mf.with_solvent.max_memory = min(mf.with_solvent.max_memory,8192)
                    mf.with_solvent.build()
                    densities = np.stack([dm,.9*dm])
                    potential = mf.with_solvent._get_v(densities)
                    q = np.random.default_rng(20261024).normal(0,.01,(2,potential.shape[1]))
                    matrix = mf.with_solvent._get_vmat(q)
                    outputs.append((potential,matrix))
                    np.savez(out/f'mult{mult}_{mode}.npz',potential=potential,matrix=matrix)
            stats = json.loads((out/f'mult{mult}_bounded/status.json').read_text())
            errors = dict(potential=float(np.max(abs(outputs[0][0]-outputs[1][0]))),matrix=float(np.max(abs(outputs[0][1]-outputs[1][1]))))
            if (any(not np.isfinite(x) or x>1e-10 for x in errors.values())
                    or stats['maximum_tile_points'] != 400 or min(stats['tiles'].values()) < 2):
                raise ValueError('The forced multi-tile contractions did not pass')
            report['cases'].append(dict(mult=mult,errors=errors,tiles=stats['tiles'],maximum_tile_points=stats['maximum_tile_points']))
            report['source_sha256'][str(checkpoint)] = api.digest(checkpoint)
        api.verify_sources(report['source_sha256'])
        report.update(status='two_spin_CPU_PCM_multitile_contractions_passed',physical_model_changed=False,
            large_case_qualified=False,output_sha256={str(p):api.digest(p) for p in out.rglob('*') if p.is_file() and p != out/'status.json'})
    except BaseException as error:
        report.update(status='CPU_PCM_controls_failed_preserved',error=repr(error))
        raise
    finally:
        report['finished']=time.time()
        atomic_write_json(out/'status.json',report)


def prepare(source, control, out):
    original = api.read_plan(source)
    api.verify_sources(original['producer_sha256'])
    api.verify_sources(original['input_sha256'])
    proof = json.loads((control/'status.json').read_text())
    if proof['status'] != 'two_spin_CPU_PCM_multitile_contractions_passed' or len(proof['cases']) != 2:
        raise ValueError('Both forced multi-tile controls must pass')
    api.verify_sources(proof['source_sha256'])
    api.verify_sources(proof['output_sha256'])
    failed = json.loads((source/'electronic/status.json').read_text())
    owner = api.load_module('CPU_PCM_original_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if failed['status'] != 'split_resource_point_failed_preserved' or owner.owner_alive(failed):
        raise ValueError('Preserve the original failure and verify its exact owner ended')
    out.mkdir(parents=True,exist_ok=False)
    inputs = dict(original['input_sha256'])
    for p in (source/'plan.json',source/'electronic/status.json',source/'electronic/resource_stop.json',control/'status.json'):
        inputs[str(p)] = api.digest(p)
    sources = dict(original['producer_sha256'],**proof['source_sha256'])
    sources.update({str(p):api.digest(p) for p in (Path(__file__).resolve(),ADAPTER)})
    plan = dict(original,input_sha256=inputs,producer_sha256=sources,CPU_PCM_memory_MiB=8192,
        CPU_PCM_controls=str(control),prior_failed_point=str(source/'electronic/status.json'),
        geometry_qualification_origin=str(source/'geometry/status.json'),prepared=time.time())
    atomic_write_json(out/'plan.json',plan)
    print(json.dumps(dict(status='bounded_PCM_large_CPU_reference_point_prepared',PCM_memory_MiB=8192)))


def run(campaign, scratch, threads):
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256'])
    worker = api.load_module('unchanged_CPU_resource_reference_worker',ROOT/'scripts/validation/qualify_split_polymer.py')
    adapter = api.load_module('qualified_CPU_PCM_budget',ADAPTER)
    with adapter.bounded_pcm_memory(campaign/'PCM_budget',plan['CPU_PCM_memory_MiB']):
        worker.run(campaign,'electronic',threads,scratch)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('controls','prepare','run'))
    for name in ('source','controls','campaign','scratch'):
        parser.add_argument('--'+name,type=Path)
    parser.add_argument('--threads',type=int,default=16)
    args=parser.parse_args()
    if args.campaign is None:parser.error('--campaign required')
    if args.action=='controls':
        if args.source is None:parser.error('--source required')
        controls(args.source.resolve(),args.campaign.resolve())
    elif args.action=='prepare':
        if args.source is None or args.controls is None:parser.error('--source and --controls required')
        prepare(args.source.resolve(),args.controls.resolve(),args.campaign.resolve())
    else:
        if args.scratch is None:parser.error('--scratch required')
        run(args.campaign.resolve(),args.scratch.resolve(),args.threads)


if __name__=='__main__':main()
