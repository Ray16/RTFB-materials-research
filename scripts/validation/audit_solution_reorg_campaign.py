#!/usr/bin/env python
"""Independent native-SMD cohort bookkeeping, density and benchmark audit."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT,DFT,HARTREE_EV
from redox.core.protocol import geom_sha1
from redox.qm.dft import _ecp_for

OUT=ROOT/'results/validation/physical_reorg_2026-10-11'


def main():
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=4);a=p.parse_args()
    queue=json.loads((OUT/'native_smd_queue.json').read_text())
    records={};checks=[];errors=[];overlaps={}
    def check(label,passed,detail=None):
        checks.append(dict(label=label,passed=bool(passed),detail=detail))
        if not passed:errors.append(label)
    for path in sorted((OUT/'native_smd_scf').glob('*.json')):
        r=json.loads(path.read_text());k=r['key'];identity=k['identity'];geom=k['geometry']
        tag=(identity['id'],identity['state'],geom['state'],
             k['reference']['state'] if k['reference'] else None)
        digest=hashlib.sha256(json.dumps(k,sort_keys=True).encode()).hexdigest()[:20]
        check(path.name+':key_address',path.stem.endswith(digest))
        check(path.name+':converged_finite',r['converged'] is True and np.isfinite([r['e_Ha'],r['e_eV']]).all())
        check(path.name+':units',abs(r['e_eV']-r['e_Ha']*HARTREE_EV)<1e-8)
        energy_guesses=[s['e_Ha'] for s in r['scf_guesses'] if s['converged']]
        check(path.name+':lowest_guess',bool(energy_guesses) and abs(r['e_Ha']-min(energy_guesses))<1e-10)
        xyz=DFT/geom['id']/geom['state']/'opt.xyz'
        check(path.name+':geometry',geom_sha1(xyz)==geom['geom_sha1'])
        for source,h in k['source'].items():
            check(path.name+':source:'+source,hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==h)
        slot=(str(xyz),k['method']['basis'])
        if slot not in overlaps:
            atoms=read(str(xyz));symbols=atoms.get_chemical_symbols()
            mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),basis=k['method']['basis'],
                ecp=_ecp_for(symbols,k['method']['basis']),charge=identity['charge'],spin=identity['mult']-1,verbose=0)
            overlaps[slot]=mol.intor('int1e_ovlp')
        with np.load(path.with_suffix('.npz')) as npz:dm=npz['dm']
        total=dm.sum(axis=0) if dm.ndim==3 else dm
        trace=float(np.einsum('ij,ji',total,overlaps[slot]))
        check(path.name+':density_trace',abs(trace-r['expected_electrons'])<1e-6,trace)
        check(path.name+':density_hermitian',np.isfinite(dm).all() and np.allclose(total,total.T,atol=1e-10))
        if tag in records:
            check(path.name+':duplicate',abs(records[tag]['e_eV']-r['e_eV'])<1e-6)
        records[tag]=r
    cases=[]
    for c in queue['cases'][:a.limit]:
        gid,o,red=c['id'],c['oxidized'],c['reduced']
        file=OUT/f'native_smd_{gid}_{o}_{red}.json'
        if not file.exists():
            cases.append(dict(id=gid,couple=f'{o}->{red}',status='pending'));continue
        r=json.loads(file.read_text());check(file.name+':case_identity',r['case']==c)
        reconstructed=[]
        for d in r['directional']:
            i,f=d['initial'],d['final']
            ieq=records[(gid,i,i,None)];feq=records[(gid,f,f,None)]
            cross=records[(gid,f,i,None)];neq=records[(gid,f,i,i)]
            calculated=dict(total_solution_relaxation_eV=neq['e_eV']-feq['e_eV'],
                equilibrium_geometry_component_eV=cross['e_eV']-feq['e_eV'],
                frozen_slow_increment_eV=neq['e_eV']-cross['e_eV'],
                vertical_charge_change_eV=neq['e_eV']-ieq['e_eV'],
                cds_relaxation_eV=(neq['cds_Ha']-feq['cds_Ha'])*HARTREE_EV)
            for quantity,value in calculated.items():check(file.name+':'+i+':'+quantity,abs(d[quantity]-value)<1e-8)
            reconstructed.append(calculated['total_solution_relaxation_eV'])
        mean=sum(reconstructed)/2
        check(file.name+':symmetric_mean',abs(mean-r['summary']['symmetric_single_body_approximation_eV'])<1e-8)
        check(file.name+':self_exchange',abs(sum(reconstructed)-r['summary']['self_exchange_infinite_separation_eV'])<1e-8)
        check(file.name+':CDS_directional_cancellation',abs(sum(d['cds_relaxation_eV'] for d in r['directional']))<1e-8)
        item=dict(id=gid,couple=f'{o}->{red}',status='completed',mean_eV=mean,
                  physical_flags=r['summary']['physical_flags'],
                  endpoint_evidence=r['endpoint_evidence'])
        if gid=='bm_2_6_dimethyl_1_4_benzoquinone':
            item['experimental_interval_eV']=[1.26,1.28]
            item['difference_to_midpoint_eV']=mean-1.27
            item['within_reported_interval']=1.26<=mean<=1.28
            item['experimental_caveat']='CVfit-dependent; no double-layer correction; no fitted model correction'
        cases.append(item)
    result=dict(checks=len(checks),failed_checks=len(errors),errors=errors,scf_records=len(records),
                expected_cases=min(a.limit,len(queue['cases'])),cases=cases,
                completed_cases=sum(c['status']=='completed' for c in cases),
                production_promotion=False)
    (OUT/'native_smd_independent_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    (OUT/'native_smd_independent_checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False))
    print(json.dumps(result))
    if errors:raise SystemExit(2)


if __name__=='__main__':main()
