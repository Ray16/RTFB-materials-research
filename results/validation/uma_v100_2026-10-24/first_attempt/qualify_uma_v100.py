#!/usr/bin/env python
"""Test an updated V100 CUDA12 stack in explicit uncompiled FP32 general mode."""
import argparse
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();verify_gate_owner()
    import torch
    if not torch.cuda.is_available() or torch.cuda.get_device_capability()!=(7,0) or 'sm_70' not in torch.cuda.get_arch_list():
        raise RuntimeError('This compatibility test requires a working Volta/V100 CUDA build')
    from fairchem import core
    from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
    original_predict=core.pretrained_mlip.get_predict_unit
    original_calculator=core.FAIRChemCalculator
    observed=[]

    def predictor(*args,**kwargs):
        kwargs['inference_settings']=InferenceSettings(nitro=False,tf32=False,merge_mole=False,
            compile=False,activation_checkpointing=True,internal_graph_gen_version=2,execution_mode='general')
        return original_predict(*args,**kwargs)

    class CheckedFP32(original_calculator):
        def __init__(self,predictor,*args,**kwargs):
            self._checked_predictor=predictor
            super().__init__(predictor,*args,**kwargs)
        def calculate(self,*args,**kwargs):
            result=super().calculate(*args,**kwargs)
            settings=self._checked_predictor.inference_settings
            if (settings.nitro or settings.tf32 or settings.compile or settings.merge_mole
                    or str(settings.execution_mode)!='general'):
                raise RuntimeError('V100 reference left the declared FP32 general mode')
            observed.append(dict(nitro=False,tf32=False,compile=False,merge_mole=False,execution_mode='general'))
            return result
    core.pretrained_mlip.get_predict_unit=predictor
    core.FAIRChemCalculator=CheckedFP32
    source=ROOT/'scripts/validation/qualify_uma_nitro.py'
    module=api.load_module('unchanged_UMA_exact_point_worker',source)
    try:
        module.run(a.campaign.resolve(),'reference',a.out.resolve())
        report=json.loads((a.out/'status.json').read_text())
        report.update(execution_profile='V100_CUDA12_general_FP32',actual_Nitro_executed=False,
            post_evaluation_mode_checks=len(observed),checked_mode=observed[-1] if observed else None,
            compatibility_wrapper_sha256=api.digest(Path(__file__).resolve()),checked_at=time.time())
        atomic_write_json(a.out/'status.json',report)
    finally:
        core.pretrained_mlip.get_predict_unit=original_predict
        core.FAIRChemCalculator=original_calculator


if __name__=='__main__':main()
