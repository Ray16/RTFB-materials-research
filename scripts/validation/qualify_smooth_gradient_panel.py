#!/usr/bin/env python
"""All-axis CPU gradient checks for imide and amine CDS environments."""
import argparse
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(1))
    from pyscf.solvent import smd_experiment
    producer=ROOT/'scripts/validation/qualify_smooth_cds.py'
    helper=api.load_module('smooth_CDS_exact_component_functions',producer)
    source=args.source.resolve();original=json.loads(source.read_text());api.verify_sources(original['source_sha256'])
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    gradient,gradient_source=helper.gradient_functions(smd_experiment,out/'experimental_gradient_source.py')
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),producer,source,Path(smd_experiment.__file__),gradient_source)}
    rows={};report=dict(status='smooth_amine_imide_full_gradient_checks_running',started=time.time(),
        states=rows,source_sha256=sources,acceptance_Ha_Bohr=1e-6,grid=2030,
        quantum_jobs_started=0,primary_records_changed=False)
    atomic_write_json(out/'summary.json',report)
    try:
        for key in ('precursor','base','protonated','pair'):
            state=Path(original['reagent_states'][key]['state_dir']);xyz=state/'opt.xyz'
            raw=json.loads((state/'result.json').read_text());sources[str(xyz)]=api.digest(xyz)
            obj=helper.solvent_object(xyz,raw['charge'],raw['mult'],2030)
            selectors=[(i,a) for i in range(obj.mol.natm) for a in range(3)]
            result=helper.sampled_gradient_check(obj,smd_experiment.get_cds,gradient,selectors)
            result.update(full_Cartesian_gradient_checked=True,state_dir=str(state),charge=raw['charge'],mult=raw['mult'])
            rows[key]=result;atomic_write_json(out/'summary.json',report)
            if result['maximum_error_Ha_Bohr']>1e-6:raise ValueError('Full component gradient failed: '+key)
        api.verify_sources(sources)
        report['status']='smooth_amine_imide_full_gradients_checked'
    except BaseException as error:
        report.update(status='smooth_full_gradient_check_failed_preserved',error=repr(error));raise
    finally:
        report['finished']=time.time();atomic_write_json(out/'summary.json',report)
    print(json.dumps(dict(status=report['status'],maximum_error_Ha_Bohr=max(v['maximum_error_Ha_Bohr'] for v in rows.values()))))


if __name__=='__main__':main()
