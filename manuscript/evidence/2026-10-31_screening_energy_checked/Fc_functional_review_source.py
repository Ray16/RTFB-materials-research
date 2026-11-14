#!/usr/bin/env python
"""Review frozen-geometry Fc energy gaps without fitting/changing a reference."""
import argparse
import json
import math
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def metrics(errors):
    return dict(bias_V=math.fsum(errors)/len(errors),MAE_V=math.fsum(abs(v) for v in errors)/len(errors),
        RMSE_V=math.sqrt(math.fsum(v*v for v in errors)/len(errors)))


def review(root):
    owner=api.load_module('Fc_diagnostic_ended_owners',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    runtime_path=ROOT/'environments/cutensor22/runtime.json';runtime_profile=json.loads(runtime_path.read_text())
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),root/'campaign.json',runtime_path)}
    methods={}
    for method in ('wb97mv','b3lyp','pbe0'):
        values={}
        for state in ('neu','ox'):
            case=root/(method+'__'+state);plan=api.read_plan(case);path=case/'electronic/status.json'
            row=json.loads(path.read_text())
            if row['status']!='frozen_geometry_Fc_electronic_control_complete' or owner.owner_alive(row['owner']):
                raise ValueError('Every reference electronic control must finish normally')
            api.verify_sources(plan['input_sha256']);api.verify_sources(plan['producer_sha256']);api.verify_sources(row['output_sha256'])
            raw_path=Path(row['raw_root'])/'electronic_diagnostic.json';raw=json.loads(raw_path.read_text())
            if (not raw['converged_smd'] or not raw['converged_gas'] or raw['geometry_reoptimized']
                    or raw['thermochemistry_computed'] or raw['potential_computed']
                    or raw['sp_xc']!=plan['task']['xc'] or raw['sp_basis']!='def2-tzvpd'
                    or raw['charge']!=plan['task']['charge'] or raw['mult']!=plan['task']['mult']):
                raise ValueError('Control protocol or diagnostic scope differs')
            if any(not math.isfinite(raw[k]) for k in ('e_smd_eV','e_gas_eV','smd_spin_contam','gas_spin_contam')):
                raise ValueError('Nonfinite reference electronic/spin data')
            if max(abs(raw[k]) for k in ('smd_spin_contam','gas_spin_contam'))>.1:
                raise ValueError('Reference spin contamination failed')
            if method=='wb97mv' and abs(raw['e_smd_eV']-plan['primary_E_smd_eV'])>1e-6:
                raise ValueError('Baseline electronic reference did not reproduce')
            rt_path=case/'tensor_runtime/status.json';rt=json.loads(rt_path.read_text())
            if (rt['status']!='task_returned_under_checked_cutensor_runtime' or rt['pid']!=row['owner']['pid']
                    or rt['versions']!=runtime_profile['versions'] or rt['profile_sha256']!=api.digest(runtime_path)
                    or rt['library_sha256']!=runtime_profile['library_sha256']):
                raise ValueError('Actual completed runtime differs from its definition')
            guard=case/'residency_guard';g=json.loads((guard/'state.json').read_text())
            if g['foreign_resident_PIDs'] or not g['own_residency_observed'] or (guard/'interruption.json').exists():
                raise ValueError('Reference diagnostic did not retain exclusive residency')
            for p in (case/'plan.json',path,raw_path,rt_path,guard/'owner.json',guard/'state.json'):proof[str(p)]=api.digest(p)
            proof.update(plan['input_sha256']);proof.update(plan['producer_sha256']);proof.update(row['output_sha256'])
            values[state]=raw
        gas=values['ox']['e_gas_eV']-values['neu']['e_gas_eV']
        smd=values['ox']['e_smd_eV']-values['neu']['e_smd_eV']
        methods[method]=dict(gas_electronic_gap_eV=gas,SMD_electronic_gap_eV=smd,
            differential_solvation_eV=smd-gas,thermochemistry_recomputed=False,full_alternate_potential_qualified=False)
    base=methods['wb97mv'];benchmark_path=ROOT/'results/validation/publication_benchmark/current.json'
    benchmark=json.loads(benchmark_path.read_text());proof[str(benchmark_path)]=api.digest(benchmark_path)
    for method,row in methods.items():
        shift=row['SMD_electronic_gap_eV']-base['SMD_electronic_gap_eV']
        row.update(SMD_reference_electronic_shift_eV=shift,
            gas_gap_shift_eV=row['gas_electronic_gap_eV']-base['gas_electronic_gap_eV'],
            solvation_gap_shift_eV=row['differential_solvation_eV']-base['differential_solvation_eV'])
        if abs(row['gas_gap_shift_eV']+row['solvation_gap_shift_eV']-shift)>1e-8:
            raise ValueError('Reference electronic decomposition failed to close')
        row['reference_only_counterfactual']=metrics([r['error_V']-shift for r in benchmark['outliers']])
    result=dict(status='frozen_geometry_Fc_functional_sensitivity_independently_checked',checked=time.time(),methods=methods,
        source_sha256=proof,source_geometries='Same qualified R2SCAN-D4/SVP/SMD neutral and cation geometries for every functional',
        basis='uniform def2-TZVPD',primary_reference_modified=False,primary_benchmark_modified=False,
        nuclear_thermochemistry_borrowed=False,alternative_potentials_calculated=False,
        interpretation='Electronic-gap sensitivity only. The reference-only counterfactual holds every primary molecular/nuclear term fixed and is not a recalibrated or validated alternative protocol. Published other-functional reference values are not method matched by functional name alone.')
    atomic_write_json(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],methods=methods)))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();review(args.root.resolve())


if __name__=='__main__':main()
