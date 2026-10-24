#!/usr/bin/env python
"""Exercise the current dispatcher or a matched-device official FP32/Nitro pair."""
import argparse
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner


def evaluate(campaign,mode,out):
    from redox.qm import uma
    from redox.qm.uma_current import calculator_for_profile
    from redox.qm.uma_model_identity import model_identity
    original_factory,original_identity=uma.make_calculator,uma.model_identity
    watched={str(p):api.digest(p) for p in (Path(__file__).resolve(),ROOT/'src/redox/qm/uma_current.py',ROOT/'config/uma_execution.json')}
    observed=[]
    restore=[]
    if mode=='legacy':
        uma.make_calculator=lambda model,device:calculator_for_profile(model,device,'legacy',legacy_factory=original_factory)
    else:
        uma.model_identity=model_identity
    if mode=='reference':
        from fairchem import core
        from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
        previous_predictor=core.pretrained_mlip.get_predict_unit
        previous_calculator=core.FAIRChemCalculator
        def predictor(*args,**kwargs):
            kwargs['inference_settings']=InferenceSettings(nitro=False,tf32=False,merge_mole=False,
                compile=False,activation_checkpointing=True,internal_graph_gen_version=2,execution_mode='general')
            return previous_predictor(*args,**kwargs)
        class Checked(previous_calculator):
            def __init__(self,predictor,*args,**kwargs):
                self._checked=predictor;super().__init__(predictor,*args,**kwargs)
            def calculate(self,*args,**kwargs):
                result=super().calculate(*args,**kwargs)
                s=self._checked.inference_settings
                if s.nitro or s.tf32 or s.merge_mole or s.compile or str(s.execution_mode)!='general':
                    raise ValueError('The reference left the declared standard FP32 general mode')
                observed.append(dict(nitro=False,tf32=False,merge_mole=False,compile=False,execution_mode='general'))
                return result
        core.pretrained_mlip.get_predict_unit=predictor;core.FAIRChemCalculator=Checked
        restore=[core,previous_predictor,previous_calculator]
    try:
        worker=api.load_module('unchanged_frozen_UMA_panel',ROOT/'scripts/validation/qualify_uma_nitro.py')
        worker.run(campaign,mode,out)
        path=out/'status.json';record=json.loads(path.read_text())
        record.update(current_dispatch_validation=(mode=='legacy'),
            cold_model_setup_timing_valid=(mode!='legacy'),
            explicit_reference_mode_checks=len(observed),
            execution_source_sha256=watched)
        api.verify_sources(watched)
        atomic_write_json(path,record)
    finally:
        uma.make_calculator,uma.model_identity=original_factory,original_identity
        if restore:
            core,previous_predictor,previous_calculator=restore
            core.pretrained_mlip.get_predict_unit=previous_predictor;core.FAIRChemCalculator=previous_calculator


def pair(campaign,out):
    verify_gate_owner()
    out.mkdir(parents=True,exist_ok=False)
    status=dict(status='matched_device_UMA_pair_running',started=time.time(),
                quantum_energy_protocol_changed=False,legacy_equivalence_assumed=False)
    atomic_write_json(out/'status.json',status)
    try:
        for mode in ('reference','nitro'):evaluate(campaign,mode,out/mode)
        from redox.validation.uma_parity import compare_points
        report=compare_points(api.read_plan(campaign),json.loads((out/'reference/status.json').read_text()),
                              json.loads((out/'nitro/status.json').read_text()))
        report['source_sha256']={str(p):api.digest(p) for p in (campaign/'plan.json',out/'reference/status.json',out/'nitro/status.json',Path(__file__).resolve())}
        atomic_write_json(out/'comparison.json',report)
        status.update(status='matched_device_UMA_pair_finished_independent_review_pending',
                      comparison=str(out/'comparison.json'))
    except BaseException as error:
        status.update(status='matched_device_UMA_pair_failed_preserved',error=repr(error));raise
    finally:
        status['finished']=time.time();atomic_write_json(out/'status.json',status)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('legacy','pair'))
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();verify_gate_owner()
    if a.action=='legacy':
        import torch
        from redox.qm.uma_current import select_backend
        if select_backend(torch.cuda.get_device_capability())!='legacy':
            raise ValueError('The default-dispatch fallback integration test requires V100')
        evaluate(a.campaign.resolve(),'legacy',a.out.resolve())
    else:pair(a.campaign.resolve(),a.out.resolve())


if __name__=='__main__':main()
