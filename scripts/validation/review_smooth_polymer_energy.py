#!/usr/bin/env python
"""Independent alternate-energy/geometry/resource join; no primary admission."""
import argparse
import json
import math
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def review(campaign, out):
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    path = campaign / 'experimental_energy/status.json'; row = json.loads(path.read_text())
    owner = api.load_module('smooth_energy_review_owner', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    if (row['status'] != 'experimental_hybrid_energy_finished_independent_review_pending'
            or owner.owner_alive(row['owner']) or row['native_G_computed'] or row['primary_candidate_admission']):
        raise ValueError('A finished alternate electronic calculation with stopped ownership is required')
    api.verify_sources(row['output_sha256'])
    raw = Path(row['raw_root']); energy_path = raw / 'experimental_energy.json'
    energy = json.loads(energy_path.read_text())
    worker = ROOT / 'scripts/pipeline/smooth_polymer_energy.py'
    api.load_module('smooth_energy_review_gates', worker).check_energy(
        energy, plan['tasks'][0], plan['experimental_method_id'])
    from redox.core.protocol import geom_sha1
    xyz = raw / 'geometry.xyz'
    if (api.digest(xyz) != plan['scored_geometry_sha256']
            or energy['energy_geom_sha1'] != geom_sha1(xyz)
            or energy['experimental_model'] != plan['experimental_model']):
        raise ValueError('Electronic energy belongs to another geometry or solvent model')
    for medium in ('smd', 'gas'):
        value = energy[f'e_{medium}_Ha'] * 27.211386245988
        if not math.isfinite(value) or abs(value-energy[f'e_{medium}_eV']) > 1e-6:
            raise ValueError('Stored electronic energy units do not reconstruct')
    if (abs(row['E_smd_eV'] - energy['e_smd_eV']) > 1e-8
            or abs(row['E_gas_eV'] - energy['e_gas_eV']) > 1e-8
            or abs(energy['e_smd_eV'] - energy['e_gas_eV'] - energy['dG_solv_eV']) > 1e-8):
        raise ValueError('Experimental electronic record arithmetic does not close')
    model_path = campaign / 'model/model.json'; model = json.loads(model_path.read_text())
    if (model['status'] != 'experimental_smooth_CDS_context_finished_QC_separate'
            or model['energy_calls'] < 1
            or any(model['profile'][k] != v for k, v in plan['experimental_model'].items())):
        raise ValueError('The declared alternate surface must be exercised by the actual energy calculation')
    api.verify_sources(model['profile']['implementation_sha256'])
    rt_path = campaign / 'tensor_runtime/status.json'; rt = json.loads(rt_path.read_text())
    profile_path = Path(plan['runtime_profile']); profile = json.loads(profile_path.read_text())
    if (rt['status'] != 'task_returned_under_checked_cutensor_runtime' or rt['pid'] != row['owner']['pid']
            or rt['profile_sha256'] != api.digest(profile_path) or rt['versions'] != profile['versions']
            or rt['library_sha256'] != profile['library_sha256'] or rt['total_GPU_bytes'] < 31*2**30):
        raise ValueError('Actual completed electronic runtime differs from its definition')
    guard = campaign / 'residency_guard'
    identity = json.loads((guard / 'owner.json').read_text())
    residency = json.loads((guard / 'state.json').read_text())
    host_path = campaign / 'host_guard/status.json'; host = json.loads(host_path.read_text())
    if (not identity['guard_stopped'] or residency['foreign_resident_PIDs']
            or not residency['own_residency_observed'] or (guard / 'interruption.json').exists()
            or host['status'] != 'program_returned_review_required'
            or not 0 < host['peak_RSS_bytes'] <= plan['maximum_host_RSS_GiB']*2**30
            or not 0 < row['observed_peak_GPU_MiB'] <= plan['maximum_peak_GPU_MiB']):
        raise ValueError('Completed electronic ownership or resource bounds failed')
    geometry_review = Path(plan['geometry_review'])
    geometry = json.loads(geometry_review.read_text()); api.verify_sources(geometry['source_sha256'])
    if geometry['portable_method_id'] != plan['experimental_method_id']:
        raise ValueError('Geometry and electronic model identities differ')
    proof = {**plan['producer_sha256'], **plan['input_sha256'], **row['output_sha256']}
    for p in (Path(__file__).resolve(), campaign / 'plan.json', path, xyz, energy_path, model_path,
              rt_path, host_path, guard / 'owner.json', guard / 'state.json', geometry_review):
        proof[str(p)] = api.digest(p)
    proof.update({str(p):api.digest(p) for p in (raw / 'SCF_checkpoints').rglob('*.json')})
    result = dict(status='alternate_electronic_geometry_and_resource_checked_native_G_pending',
        checked=time.time(), campaign=str(campaign), task=plan['tasks'][0],
        experimental_model=plan['experimental_model'], portable_method_id=plan['experimental_method_id'],
        E_smd_eV=energy['e_smd_eV'], E_gas_eV=energy['e_gas_eV'],
        observed_peak_GPU_MiB=row['observed_peak_GPU_MiB'], peak_host_RSS_bytes=host['peak_RSS_bytes'],
        source_sha256=proof, native_Hessian_or_G_qualified=False, primary_candidate_admission=False,
        scope='Independently stationary alternate-model geometry and fresh hybrid energy; native thermochemistry/ensemble separate')
    out.parent.mkdir(parents=True, exist_ok=True); atomic_write_json(out, result)
    print(json.dumps({k:result[k] for k in ('status', 'E_smd_eV', 'observed_peak_GPU_MiB')}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(); review(args.campaign.resolve(), args.out.resolve())


if __name__ == '__main__': main()
