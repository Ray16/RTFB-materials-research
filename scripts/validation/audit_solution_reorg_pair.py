#!/usr/bin/env python
"""Independent density/directional-energy audit; optional historical AQ comparison."""
import hashlib
import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT,HARTREE_EV
from redox.qm.dft import _ecp_for
from redox.core.protocol import load_record

BASE=ROOT/'results/validation/physical_reorg_2026-10-11'
OUT=BASE/'graft_local_solution_pair'


def main():
    global OUT
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=OUT)
    parser.add_argument('--report-prefix',default='graft_local_pair_independent')
    parser.add_argument('--report-dir',type=Path,default=BASE)
    parser.add_argument('--pair-only',action='store_true',help='audit a self-contained pair without the historical anthraquinone comparison')
    args=parser.parse_args();OUT=args.out
    result_bytes=(OUT/'result.json').read_bytes()
    result=json.loads(result_bytes);identity=result['identity'];states=identity['pair']['states']
    evidence_paths=set(OUT.glob('scf_*.json'))|set(OUT.glob('scf_*.npz'))
    evidence_paths.update(Path(s['geometry']) for s in states)
    evidence_paths.update(Path(s) for s in identity['source_sha256'])
    evidence_paths.update(Path(s) for s in identity['pair'].get('implementation_controls',{}))
    evidence_paths.update(Path(e['minimum_source']) for e in result['endpoint_evidence'] if e.get('minimum_source'))
    evidence={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in evidence_paths}
    evidence[str((OUT/'result.json').resolve())]=hashlib.sha256(result_bytes).hexdigest()
    checks=[]
    def check(label,ok,detail=None):checks.append(dict(label=label,passed=bool(ok),detail=detail))
    by_geometry={}
    for s in states:
        data=Path(s['geometry']).read_bytes();sha256=hashlib.sha256(data).hexdigest()
        check(s['name']+':source_geometry',sha256==identity['input_sha256'][s['geometry']])
        by_geometry[sha256]=s
    current={};paths={}
    for path in OUT.glob('scf_*.json'):
        r=json.loads(path.read_text());key=r['key'];state=key['state'];geometry=by_geometry[key['geometry_sha256']]
        digest=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest()[:20]
        check(path.name+':key_address',path.stem=='scf_'+digest)
        check(path.name+':pair_identity',key['identity']==identity)
        check(path.name+':converged_finite',r['converged'] is True and np.isfinite([r['e_Ha'],r['e_eV']]).all())
        atoms=read(geometry['geometry']);symbols=atoms.get_chemical_symbols()
        mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),charge=state['charge'],spin=state['mult']-1,
                  basis=identity['method']['basis'],ecp=_ecp_for(symbols,identity['method']['basis']),verbose=0)
        with np.load(path.with_suffix('.npz')) as npz:dm=npz['dm']
        total=dm.sum(axis=0) if dm.ndim==3 else dm
        trace=float(np.einsum('ij,ji',total,mol.intor('int1e_ovlp')))
        check(path.name+':electron_trace',abs(trace-mol.nelectron)<1e-6,trace)
        check(path.name+':finite_density',np.isfinite(dm).all())
        check(path.name+':hermitian_density',np.allclose(total,total.T,atol=1e-10))
        check(path.name+':energy_units',abs(r['e_eV']-r['e_Ha']*HARTREE_EV)<1e-8)
        check(path.name+':lowest_converged',abs(r['e_Ha']-min(g['e_Ha'] for g in r['scf_guesses'] if g['converged']))<1e-10)
        tag=(state['name'],geometry['name'],key['reference'] is not None)
        check(path.name+':unique_leg',tag not in current)
        current[tag]=r;paths[tag]=path
    for tag,r in current.items():
        reference=r['key']['reference']
        if reference is None:continue
        ref_tag=(reference['name'],reference['name'],False)
        ref_path=paths.get(ref_tag)
        check(str(tag)+':reference_state',reference in states and tag[1]==reference['name'])
        if ref_path is not None:
            with np.load(ref_path.with_suffix('.npz')) as npz:reference_hash=hashlib.sha256(npz['dm'].tobytes()).hexdigest()
            check(str(tag)+':reference_density',r['key']['reference_density_sha256']==reference_hash)
        elif args.pair_only:check(str(tag)+':reference_density_present',False)
    neutral=states[0];red=states[1]
    neutral_sha1=hashlib.sha1(Path(neutral['geometry']).read_bytes()).hexdigest()[:16]
    inherited={}
    for path in (() if args.pair_only else (BASE/'native_smd_scf').glob('*.json')):
        r=json.loads(path.read_text());k=r['key']
        if k['geometry']['geom_sha1']!=neutral_sha1 or k['method']!=identity['method']:continue
        q=k['identity']['charge'];ref=k['reference'] is not None
        inherited[(q,ref)]=r
    def selected(species,geometry,reference):
        key=(species['name'],geometry['name'],reference)
        if key in current:return current[key]
        if geometry['name']!=neutral['name']:raise RuntimeError('No audited SCF for requested leg')
        return inherited[(species['charge'],reference)]
    n=selected(neutral,neutral,False);r_cross=selected(red,neutral,False)
    r_neq=selected(red,neutral,True);r_eq=selected(red,red,False)
    n_cross=selected(neutral,red,False);n_neq=selected(neutral,red,True)
    expected=[dict(total_solution_relaxation_eV=r_neq['e_eV']-r_eq['e_eV'],
        equilibrium_geometry_component_eV=r_cross['e_eV']-r_eq['e_eV'],
        frozen_slow_increment_eV=r_neq['e_eV']-r_cross['e_eV'],
        vertical_charge_change_eV=r_neq['e_eV']-n['e_eV']),
        dict(total_solution_relaxation_eV=n_neq['e_eV']-n['e_eV'],
        equilibrium_geometry_component_eV=n_cross['e_eV']-n['e_eV'],
        frozen_slow_increment_eV=n_neq['e_eV']-n_cross['e_eV'],
        vertical_charge_change_eV=n_neq['e_eV']-r_eq['e_eV'])]
    for i,(actual,values) in enumerate(zip(result['directional'],expected)):
        for quantity,value in values.items():check(f'direction{i}:{quantity}',abs(actual[quantity]-value)<1e-8)
    mean=sum(r['total_solution_relaxation_eV'] for r in expected)/2
    check('symmetric_average',abs(result['summary']['symmetric_single_body_approximation_eV']-mean)<1e-8)
    check('CDS_charge_independence_at_local_geometry',abs(r_eq['cds_Ha']-n_cross['cds_Ha'])<1e-12)
    check('CDS_response_independence_at_local_geometry',abs(n_neq['cds_Ha']-n_cross['cds_Ha'])<1e-12)
    neutral_active=load_record(Path(neutral['geometry']).parent,neutral['charge'],neutral['mult'])
    if neutral_active:
        check('neutral_active_production_energy_replay',abs(n['e_eV']-neutral_active['e_smd_eV'])<1e-4,
              n['e_eV']-neutral_active['e_smd_eV'])
    if args.pair_only:
        for source,digest in identity['source_sha256'].items():
            check('producer_source:'+source,hashlib.sha256(Path(source).read_bytes()).hexdigest()==digest)
        for source,digest in identity['pair'].get('implementation_controls',{}).items():
            check('implementation_control:'+source,hashlib.sha256(Path(source).read_bytes()).hexdigest()==digest)
        for source,digest in evidence.items():
            check('unchanged_during_audit:'+source,hashlib.sha256(Path(source).read_bytes()).hexdigest()==digest)
        report=dict(checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),
            fresh_SCF_records=len(current),mean_eV=mean,physical_flags=result['summary']['physical_flags'],
            evidence_sha256=evidence,raw_outputs_modified=False,production_promotion=False)
        args.report_dir.mkdir(parents=True,exist_ok=True)
        for name,data in [('audit',report),('checks',checks)]:
            (args.report_dir/f'{args.report_prefix}_{name}.json').write_text(json.dumps(data,indent=2,allow_nan=False))
        print(json.dumps(report))
        if report['failed_checks']:raise SystemExit(2)
        return
    global_case=json.loads((BASE/'native_smd_anthraquinone_neu_red1.json').read_text())
    global_mean=global_case['summary']['symmetric_single_body_approximation_eV']
    report=dict(checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),
        fresh_SCF_records=len(current),mean_eV=mean,global_mean_eV=global_mean,
        global_comparison_convention='historical original global snapshot; not a verified same-source conformer-only difference',
        local_minus_global_eV=mean-global_mean,local_minus_global_red1_electronic_eV=(r_eq['e_eV']-
            next(json.loads(p.read_text())['e_eV'] for p in (BASE/'native_smd_scf').glob('*.json')
                 if (lambda k:k['identity']['id']=='anthraquinone' and k['identity']['state']=='red1'
                     and k['geometry']['state']=='red1' and k['reference'] is None)(json.loads(p.read_text())['key']))),
        raw_outputs_modified=False,production_promotion=False)
    (BASE/f'{args.report_prefix}_audit.json').write_text(json.dumps(report,indent=2))
    (BASE/f'{args.report_prefix}_checks.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(report))
    if report['failed_checks']:raise SystemExit(2)


if __name__=='__main__':main()
