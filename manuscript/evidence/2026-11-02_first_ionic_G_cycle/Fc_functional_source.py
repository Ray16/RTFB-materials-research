#!/usr/bin/env python
"""Frozen-geometry Fc electronic-gap diagnostic; primary reference unchanged.

Rescores the same two geometries with three declared functionals. No alternate
potential or thermochemistry is computed by borrowing nuclear corrections.
"""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT,DFT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path,ACTIVE_SP
from redox.workflows import polymer_campaign as api


def prepare(out):
    if out.exists():raise ValueError('Preserve an existing sensitivity campaign')
    out.mkdir(parents=True)
    profile=ROOT/'environments/cutensor22/runtime.json'
    tasks=[]
    for method,xc,nlc in [('wb97mv','wb97m-v','vv10'),('b3lyp','b3lyp',None),('pbe0','pbe0',None)]:
        for state,charge,mult in [('neu',0,1),('ox',1,2)]:
            source=DFT/'ferrocene'/state;parent=load_record(source,charge,mult)
            if parent is None:raise ValueError('Qualified baseline Fc electronic records required')
            case=out/(method+'__'+state);case.mkdir()
            xyz=case/'input.xyz';xyz.write_bytes((source/'opt.xyz').read_bytes())
            inputs={str(p):api.digest(p) for p in (source/'opt.xyz',record_path(source,charge,mult),xyz,profile)}
            producers={str(p):api.digest(p) for p in (Path(__file__).resolve(),ROOT/'src/redox/qm/dft.py',
                ROOT/'src/redox/core/protocol.py',ROOT/'scripts/pipeline/smd_gc_workspace.py',
                ROOT/'scripts/fleet/tensor_runtime_worker.py')}
            plan=dict(task=dict(key=case.name,charge=charge,mult=mult,state=state,xc=xc,nlc=nlc,basis='def2-tzvpd'),
                source=str(source),input_xyz=str(xyz),primary_E_smd_eV=parent['e_smd_eV'],
                input_sha256=inputs,producer_sha256=producers,runtime_profile=str(profile),
                execution_driver=str(Path(__file__).resolve()),maximum_host_RSS_GiB=64,
                raw_root=str(ROOT/'calcs/validation'/out.name/case.name),
                scope='Legacy-SMD fixed-geometry electronic-only functional sensitivity; no new potential/G',
                primary_reference_modified=False,thermochemistry_borrowed=False)
            atomic_write_json(case/'plan.json',plan)
            tasks.append(dict(key=case.name,campaign=str(case),plan_sha256=api.digest(case/'plan.json')))
    atomic_write_json(out/'campaign.json',dict(status='six_electronic_reference_controls_prepared',tasks=tasks,
        primary_reference_modified=False,thermochemistry_recomputed=False))
    print(json.dumps(dict(prepared=len(tasks),primary_reference_modified=False)))


def run(campaign):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();plan=api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    task=plan['task'];raw=Path(plan['raw_root']);raw.mkdir(parents=True,exist_ok=False)
    status=campaign/'electronic';status.mkdir(exist_ok=False)
    row=dict(status='frozen_geometry_Fc_electronic_control_running',started=time.time(),task=task,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        primary_reference_modified=False,geometry_reoptimized=False,thermochemistry_computed=False,
        potential_computed=False,raw_root=str(raw))
    atomic_write_json(status/'status.json',row)
    try:
        from redox.qm import dft
        from gpu4pyscf.solvent import smd
        if smd.get_cds_legacy.__module__!='gpu4pyscf.solvent.smd':
            raise ValueError('Reference diagnostic requires the unchanged legacy CDS model')
        memory=api.load_module('Fc_reference_memory',ROOT/'scripts/pipeline/run_pi_stacking.py')
        adapter=api.load_module('Fc_reference_controlled_workspace',ROOT/'scripts/pipeline/smd_gc_workspace.py')
        with memory.MemoryMonitor() as observed,adapter.compact_smd_execution(raw/'execution',28),\
                (raw/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
            result=dft.dft_smd(Path(plan['input_xyz']),task['charge'],task['mult'],backend='gpu',
                do_opt=False,do_freq=False,do_smd=True,do_gas=True,
                sp_xc=task['xc'],sp_basis=task['basis'],sp_basis_anion=task['basis'],sp_disp=None,sp_nlc=task['nlc'])
        if not result['converged_smd'] or not result['converged_gas']:
            raise ValueError('Reference diagnostic SCFs failed')
        if not 0<observed.peak<=31744:raise ValueError('Reference diagnostic exceeded 32 GB class ceiling')
        if task['xc']=='wb97m-v' and abs(result['e_smd_eV']-plan['primary_E_smd_eV'])>1e-6:
            raise ValueError('Baseline repeat differs; resolve before interpreting other functionals')
        result.update(electronic_reference_method_control=True,primary_reference_modified=False,
            geometry_reoptimized=False,thermochemistry_computed=False,potential_computed=False)
        atomic_write_json(raw/'electronic_diagnostic.json',result)
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        row.update(status='frozen_geometry_Fc_electronic_control_complete',observed_peak_GPU_MiB=observed.peak,
            E_smd_eV=result['e_smd_eV'],E_gas_eV=result['e_gas_eV'],
            output_sha256={str(p):api.digest(p) for p in raw.rglob('*') if p.is_file() and p.suffix in ('.json','.py')})
    except BaseException as error:
        row.update(status='Fc_electronic_control_failed_preserved',error=repr(error));raise
    finally:
        row['finished']=time.time();atomic_write_json(status/'status.json',row)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','run'))
    parser.add_argument('--campaign',type=Path,required=True);args=parser.parse_args()
    if args.action=='prepare':prepare(args.campaign.resolve())
    else:run(args.campaign.resolve())


if __name__=='__main__':main()
