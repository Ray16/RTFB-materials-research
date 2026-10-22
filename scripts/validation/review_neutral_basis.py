#!/usr/bin/env python
"""Independent-input audit of the isolated neutral-SVPD panel; no QM submitted."""
import argparse,hashlib,importlib.util,json,math,os
from pathlib import Path
import numpy as np
from redox.core.common import ROOT
from redox.core.protocol import load_record,record_path
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False);plan=json.loads((args.campaign/'plan.json').read_text())
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    source={str(args.campaign/'plan.json'):sha(args.campaign/'plan.json'),str(Path(__file__)):sha(__file__)}
    spec=importlib.util.spec_from_file_location('basis_independent_partition',ROOT/'scripts/validation/audit/native_partition.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    spec=importlib.util.spec_from_file_location('basis_independent_primary',ROOT/'scripts/validation/audit/recompute_axes.py')
    primary_audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(primary_audit)
    rows=[];pending=[]
    for task in plan['tasks']:
        sd=Path(plan['raw_root'])/task['key'];status=sd/'status.json'
        if not status.exists() or json.loads(status.read_text())['status']!='alternative_basis_control_checked':
            pending.append(task['key']);continue
        data=json.loads(status.read_text());native=Path(data['surface_record']);surface=json.loads(native.read_text());geo=sd/'opt.xyz'
        if surface['identity']['protocol']['basis']!='def2-svpd' or surface['primary_native_certificate_claimed'] is not False:
            raise ValueError('Alternative-basis identity/certificate is wrong')
        gradient=surface['gradient_max'];asym=surface['maximum_asymmetry']
        if not 0<=gradient<=.00045 or not 0<=asym<=1e-6:raise ValueError('Alternative physical gate failed')
        symbols=[r.split()[0] for r in geo.read_text().splitlines()[2:]];oxygen=next((i for i,s in enumerate(symbols) if s=='O'),len(symbols)-1)
        required={(a,c,h) for a,c in ((0,1),(oxygen,2)) for h in (.0005,.001)}
        checks=surface['derivative_checks']
        if required!={(c['atom_index'],c['axis'],c['step_Bohr']) for c in checks} or any(not 0<=c['relative_error']<=.01 for c in checks):
            raise ValueError('Alternative Hessian column checks incomplete')
        matrix=np.load(native.with_suffix('')/'hessian_cartesian.npy')
        if not np.isfinite(matrix).all() or abs(float(abs(matrix-matrix.T).max())-asym)>1e-12:
            raise ValueError('Saved alternative Hessian differs')
        partition=audit.read_partition(native,geo)
        if partition is None:raise ValueError('Alternative nuclear partition failed independent reconstruction')
        source.update(partition['evidence_sha256'])
        sp=load_record(sd,0,1)
        original=Path(task['primary']['state_dir']);raw=json.loads((original/'result.json').read_text())
        original_sp=load_record(original,0,1);nuclear,nimag,definition=primary_audit.thermal_inputs(original,raw)
        if sp is None or original_sp is None or nuclear is None or nimag!=0:raise ValueError('Backing active/primary records incomplete')
        # The independent partition reader returns its nuclear/spin sum in g_thermal_total_eV.
        thermal=partition['g_thermal_total_eV'];total=sp['e_smd_eV']+thermal
        deltas=dict(delta_thermal_eV=thermal-nuclear,delta_geometry_SP_eV=sp['e_smd_eV']-original_sp['e_smd_eV'],
            delta_first_reduction_V_holding_anion_and_Fc_fixed=total-original_sp['e_smd_eV']-nuclear)
        if any(abs(deltas[k]-data[k])>1e-8 for k in deltas):raise ValueError('Reported basis comparison differs')
        for p in (status,native,geo,native.with_suffix('')/'hessian_cartesian.npy',record_path(sd,0,1),original/'result.json',record_path(original,0,1)):
            source[str(p)]=sha(p)
        rows.append(dict(id=task['key'],**deltas,gradient_max_Ha_Bohr=gradient,maximum_Hessian_asymmetry=asym,
            maximum_sampled_column_error=max(c['relative_error'] for c in checks),alternative_surface=str(native),primary_state_dir=str(original)))
    if any(sha(p)!=h for p,h in source.items()):raise ValueError('Review inputs changed')
    result=dict(status='five_neutral_basis_controls_independently_reconstructed' if len(rows)==len(plan['tasks']) else 'partial_basis_controls_checked',
        checked=len(rows),expected=len(plan['tasks']),pending=pending,rows=rows,source_sha256=source,
        primary_values_modified=False,full_redox_basis_convergence_claimed=False)
    (args.out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


if __name__=='__main__':main()
