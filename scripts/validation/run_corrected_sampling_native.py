#!/usr/bin/env python3
"""Execute the original smooth-model native driver with corrected screening.

Full native equations, sampled derivatives, mode projection, own nuclear terms
and the exact previously checked hybrid electronic input remain mandatory.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();case=a.plan.parent;sources=case.parent/'sources';plan=json.loads(a.plan.read_text())
    if (plan['role']!='corrected_sampling_native' or not plan['corrected_sampling_native']
            or plan['native_G_qualified'] or plan['primary_candidate_admission'] or plan['source_nuclear_or_G_transferred']):
        raise ValueError('Fresh corrected native method-control scope required')
    def verify():
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h for p,h in plan['source_sha256'].items()):
            raise ValueError('Frozen native source changed')
    verify();owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    stage=dict(status='corrected_sampling_native_running',owner=owner,started=time.time(),native_G_qualified=False)
    atomic_write_json(case/'status.json',stage)
    try:
        from redox.qm import dft
        from gpu4pyscf.hessian import rhf
        original=load(ROOT/'scripts/pipeline/smooth_polymer_native.py','sampling_original_native_equations')
        if original.NATIVE_SETTINGS!=plan['experimental_native_settings']:
            raise ValueError('Original native physical settings changed')
        helper=load(sources/'helper.py','sampling_native_corrected_screening')
        checkpoints=load(sources/'checkpoints.py','sampling_native_full_checkpoint')
        frames=load(sources/'native_checkpoints.py','sampling_native_grid_checkpoint')
        response=load(sources/'response_helper.py','sampling_native_response_budget')
        with helper.corrected_screening(sources/'kernel.cu',case/'screening_implementation'), \
                frames.recording_native(dft,checkpoints,Path(plan['raw_root']),owner), \
                response.response_memory_cap(rhf,4*2**30):
            original.run(case)
        verify();result=json.loads((case/'experimental_native/status.json').read_text())
        if result['status']!='experimental_native_finished_independent_review_pending':
            raise ValueError('Original native calculation did not finish normally')
        stage.update(status='corrected_sampling_native_finished_independent_review_pending',
            observed_peak_GPU_MiB=result['observed_peak_GPU_MiB'])
    except BaseException as error:
        stage.update(status='corrected_sampling_native_failed_preserved',error=repr(error));raise
    finally:
        stage['finished']=time.time();atomic_write_json(case/'status.json',stage)


if __name__=='__main__':main()
