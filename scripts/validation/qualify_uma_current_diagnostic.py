#!/usr/bin/env python
"""Replay the failed first UMA case with durable per-evaluation diagnostics.

The original scientific plan, thresholds, calculator and reference settings are
retained. Every call is recorded before the existing repeatability gate can fail.
This diagnostic does not qualify the full panel or change the default calculator.
"""
import argparse
import inspect
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ORIGINAL=ROOT/'scripts/validation/qualify_uma_current.py'
PANEL=ROOT/'scripts/validation/qualify_uma_nitro.py'


def instrument(worker,diagnostic_source):
    body=inspect.getsource(worker.run)
    old="""                return dict(energy_eV=energy,forces_eV_A=forces.tolist(),seconds=elapsed)"""
    new="""                value = dict(energy_eV=energy,forces_eV_A=forces.tolist(),seconds=elapsed)
                call = len(_calls)
                _calls.append(dict(case_id=case['id'],call=call,frame_index=index,**value))
                atomic_write_json(out/'evaluation_diagnostics.json', dict(status='per_evaluation_diagnostic_unqualified', calls=_calls, source_sha256=_diagnostic_sources))
                return value"""
    if body.count(old)!=1:raise ValueError('UMA point-evaluation template changed')
    body=body.replace(old,new,1)
    old="            if (error_E>plan['repeat_energy_tolerance_eV_per_atom']*case['atoms'] or error_F>plan['repeat_force_tolerance_eV_A']):"
    new="""            atomic_write_json(out/'repeat_diagnostic.json', dict(case_id=case['id'], mode=mode, repeat_energy_error_eV=error_E, repeat_force_error_eV_A=error_F, energy_tolerance_eV=plan['repeat_energy_tolerance_eV_per_atom']*case['atoms'], force_tolerance_eV_A=plan['repeat_force_tolerance_eV_A'], calls=len(_calls), source_sha256=_diagnostic_sources, qualification_inferred=False))
"""+old
    if body.count(old)!=1:raise ValueError('UMA repeatability gate template changed')
    body=body.replace(old,new,1)
    proof={str(p):api.digest(p) for p in (diagnostic_source,ORIGINAL,PANEL)}
    namespace=dict(vars(worker),_calls=[],_diagnostic_sources=proof);local={}
    exec(compile(body,str(diagnostic_source),'exec'),namespace,local)
    original_run=local['run']
    def run(campaign,mode,out,case_id=None):
        plan=api.read_plan(campaign)
        if case_id is not None:raise ValueError('Diagnostic uses the exact first failed case')
        api.verify_sources(proof)
        return original_run(campaign,mode,out,plan['cases'][0]['id'])
    worker.run=run
    return worker


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('pair',))
    parser.add_argument('--campaign',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    original_load=api.load_module
    def load(name,path):
        module=original_load(name,path)
        if Path(path).resolve()==PANEL:instrument(module,Path(__file__).resolve())
        return module
    current=original_load('official_UMA_first_case_diagnostic',ORIGINAL)
    api.load_module=load
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    status=dict(status='first_case_repeat_diagnostic_running',started=time.time(),
                full_panel_qualified=False,production_calculator_qualified=False)
    try:
        for mode in ('reference','nitro'):
            current.evaluate(a.campaign.resolve(),mode,out/mode)
        status['status']='first_case_repeat_diagnostic_completed_review_pending'
    except BaseException as error:
        status.update(status='first_case_repeat_diagnostic_failed_preserved',error=repr(error))
        raise
    finally:
        status['finished']=time.time();atomic_write_json(out/'status.json',status)
        api.load_module=original_load


if __name__=='__main__':main()
