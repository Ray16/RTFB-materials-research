#!/usr/bin/env python
"""Check a saved trimer geometry/SP without repeating optimization or gas SCF."""
import argparse,json,os,time
from pathlib import Path
import numpy as np
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--plan',type=Path,required=True);args=parser.parse_args()
    verify_gate_owner();plan=json.loads(args.plan.read_text());api.verify_sources(plan['source_sha256'])
    if len(os.sched_getaffinity(0))!=2:raise ValueError('Two assigned CPU slots required')
    sd=Path(plan['state_dir']);out=Path(plan['raw_out']);out.mkdir(parents=True,exist_ok=False)
    status=dict(status='checking_saved_checkpoint',started=time.time(),host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),
        optimized_geometry_recomputed=False,free_energy_computed=False,minimum_Hessian_tested=False)
    atomic_write_json(out/'status.json',status)
    try:
        from redox.validation.protonation import hydrogen_topology_status
        from redox.build.geometry import geometry_stereo_status
        from redox.build.polymer_contacts import contact_metrics
        from ase.io import read
        energy=load_record(sd,0,1)
        if energy is None:raise ValueError('Saved active SP no longer admits')
        topology=hydrogen_topology_status(plan['task']['smiles'],sd/'opt.xyz')
        stereo=geometry_stereo_status(plan['task']['smiles'],sd/'opt.xyz')
        if topology!='ok' or stereo not in ('ok','not_applicable'):raise ValueError('Saved graph/stereochemistry fails')
        returned=json.loads((sd/'geometry_checkpoints/optimizer_returned.json').read_text())
        if returned['optimizer_converged'] is not True:raise ValueError('Optimizer did not return convergence')
        adaptive=api.load_module('saved_trimer_adaptive_forces',ROOT/'scripts/pipeline/resume_polymer_forces_adaptive.py')
        with adaptive.adaptive_forces(out/'CDS_checks.jsonl'):
            force=adaptive.base.final_energy_check(sd/'opt.xyz',out/'force_energy_check')
        baseline=json.loads((out/'force_energy_check/baseline.json').read_text());gmax=float(abs(np.asarray(baseline['gradient_Ha_Bohr'])).max())
        if gmax>.00045:raise ValueError('Saved geometry is not stationary')
        atoms=read(sd/'opt.xyz');contacts=contact_metrics(atoms.positions,plan['task']['pi_groups'])
        api.verify_sources(plan['source_sha256'])
        status.update(status='saved_stationary_electronic_checkpoint_checked',E_smd_eV=energy['e_smd_eV'],
            gradient_max_Ha_Bohr=gmax,force_energy_check=force,topology_QC=topology,stereochemistry_QC=stereo,contact_metrics=contacts,
            state_dir=str(sd),geometry_sha256=api.digest(sd/'opt.xyz'),original_worker_exit=139,
            original_evidence_preserved=True,synthesis_admission=False)
    except BaseException as error:
        status.update(status='failed_or_interrupted_preserving_records',error=repr(error));raise
    finally:
        status['finished']=time.time();atomic_write_json(out/'status.json',status)


if __name__=='__main__':main()
