#!/usr/bin/env python
"""Review one completed external-control UMA resource pilot on the CPU."""
import argparse
import csv
import json
import math
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm import uma
from redox.workflows import polymer_campaign as api


def review(campaign, attempt, out):
    adapter = api.load_module('review_external_UMA_control_adapter',ROOT/'scripts/pipeline/uma_external_control.py')
    plan, rows = adapter.validate_campaign(campaign)
    row = next(r for r in rows if r['state'] == 'neu')
    directory = campaign/'output/nitro'/plan['group']/'neu'
    result_path, xyz = directory/'result.json', directory/'relaxed.xyz'
    receipt_path = campaign/'completed/neu.json'
    result, receipt = json.loads(result_path.read_text()), json.loads(receipt_path.read_text())
    if receipt['result_sha256'] != api.digest(result_path) or receipt['status'] != 'UMA_control_worker_returned_review_required':
        raise ValueError('Worker receipt is missing or differs from result')
    policy_path = ROOT/'config/uma_execution.json'
    policy = json.loads(policy_path.read_text())
    identity = result['model_identity']
    if (identity.get('checkpoint_sha256') != policy['checkpoint_sha256']
            or identity.get('reference_sha256') != policy['reference_sha256']
            or result.get('UMA_execution',{}).get('actual_mode') != 'nitro'
            or result['UMA_execution'].get('policy_sha256') != api.digest(policy_path)):
        raise ValueError('Actual model identity or Nitro execution differs')
    seeds = sorted((Path(plan['library'])/plan['group']/'conformers').glob('conf_*.xyz'))
    errors = uma.cache_errors(result,row,xyz,seeds,identity,0.02)
    if errors:
        raise ValueError('UMA physical/cache checks failed: ' + ', '.join(errors))
    guard = attempt/'residency_guard'; host_path = attempt/'host_guard/status.json'; gpu_path=attempt/'gpu_memory/status.json'
    owner = json.loads((guard/'owner.json').read_text()); state = json.loads((guard/'state.json').read_text())
    host = json.loads(host_path.read_text()); memory = json.loads(gpu_path.read_text())
    if (not owner.get('guard_stopped') or state.get('foreign_resident_PIDs') or not state.get('own_residency_observed')
            or (guard/'interruption.json').exists() or (guard/'observation_error.json').exists()
            or host.get('status') != 'program_returned_review_required'
            or not 0 < host.get('peak_RSS_bytes',0) <= 16*2**30
            or memory.get('status') != 'GPU_memory_observation_program_returned_review_required'
            or memory.get('own_samples',0) < 1 or memory.get('sampling_errors')
            or not isinstance(memory.get('peak_process_memory_MiB'),(int,float))
            or not math.isfinite(memory['peak_process_memory_MiB']) or memory['peak_process_memory_MiB'] <= 0
            or owner.get('pid') != memory.get('pid') or str(owner.get('GPU')) != str(memory.get('GPU'))
            or receipt.get('owner',{}).get('pid') != owner.get('pid')
            or receipt.get('owner',{}).get('process_start') != owner.get('process_start')
            or host.get('pid') != owner.get('pid') or host.get('CPU_affinity') != owner.get('CPU_affinity')):
        raise ValueError('Pilot residency or measured resource evidence is incomplete')
    paths = [campaign/'plan.json',receipt_path,result_path,xyz,policy_path,guard/'owner.json',guard/'state.json',
             host_path,gpu_path,Path(__file__).resolve(),ROOT/'src/redox/qm/uma.py']
    report = dict(status='NDI_neutral_UMA_pilot_independently_checked',state='neu',
        actual_backend='nitro',energy_eV=result['energy_eV'],maximum_force_eV_A=result['fmax'],
        peak_host_RSS_bytes=host['peak_RSS_bytes'],peak_process_GPU_MiB=memory['peak_process_memory_MiB'],
        GPU_memory_samples=memory['own_samples'],source_sha256={str(p):api.digest(p) for p in paths},
        solution_DFT_or_thermochemistry_qualified=False,experimental_accuracy_inferred=False,
        scope='One neutral gas-phase UMA warm-start/resource pilot; charged states and full optimizer validation remain open')
    if out.exists():
        raise FileExistsError('Preserve existing pilot review')
    atomic_write_json(out,report)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--attempt',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    result=review(a.campaign.resolve(),a.attempt.resolve(),a.out.resolve())
    print(json.dumps({'status':result['status'],'peak_process_GPU_MiB':result['peak_process_GPU_MiB']}))
