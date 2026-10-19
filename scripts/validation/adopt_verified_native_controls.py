#!/usr/bin/env python
"""Adopt independently derivative-checked controls into native thermal caches.

No optimizer JSON, geometry or electronic energy is changed. This avoids rerunning
the expensive Hessians already verified against complete numerical references.
"""
import hashlib,json,shutil,time
from pathlib import Path
import numpy as np
from redox.core.common import ROOT,DFT,read_result
from redox.qm.native_thermal import NATIVE_THERMAL,cache_path,load_native_thermal,implementation_controls
from redox.qm.thermal import spin_degeneracy_eV


def main():
    base=ROOT/'results/validation/physical_thermal_spin_2026-10-10';proof=implementation_controls();adoptions=[]
    specs=[('aq_benzylamino','red2',base/'analytic_full_response_aq_reference','analytic_full_response_verified'),
        ('bm_1_4_benzoquinone','red1',base/'uks_bq_full_response_control','uks_full_response_verified')]
    fc=base/'uks_fc_full_response_control'
    if (fc/'status.json').exists() and json.loads((fc/'status.json').read_text())['status']=='uks_full_response_verified':
        specs.append(('ferrocene','ox',fc,'uks_full_response_verified'))
    for gid,state,source,expected in specs:
        raw=read_result(gid,state,raw=True);q,m=raw['charge'],raw['mult'];sd=DFT/gid/state
        p,identity=cache_path(sd,q,m)
        if load_native_thermal(sd,q,m):continue
        if p.exists() or p.with_suffix('').exists():raise RuntimeError('preserve/inspect an existing partial native cache')
        evidence=json.loads((source/'status.json').read_text());assert evidence['status']==expected
        if evidence.get('cphf_grid_level')!=NATIVE_THERMAL.get('cphf_grid_level'):
            raise RuntimeError('Control used a different/unsaved response grid; preserve its older protocol, do not relabel it')
        gh=evidence.get('source_sha256',evidence.get('geometry_sha256'))
        assert gh==identity['geometry_sha256']
        if expected=='analytic_full_response_verified':
            values=evidence['thermal_by_cutoff_eV'];freq=evidence['frequencies_cm']
            derivative=[dict(column=r['column'],relative_error=r['relative_norm_error'],
                reference='complete conservative numerical Hessian') for r in evidence['column_checks']]
            spin_diag=dict(s_squared=0.,spin_contam=0.,source='restricted singlet identity')
            matrix=np.load(source/'hessian.npy').transpose(0,2,1,3).reshape(len(freq)+6,len(freq)+6)
            logs=source
        else:
            values=evidence['qrrho']['analytical'];freq=evidence['analytical_frequencies_cm']
            derivative=[dict(column=r['column'],relative_error=r['step_and_cold_relative_error'],
                reference='complete conservative numerical Hessian plus cold/two-step check') for r in evidence['step_cold_checks']]
            spin_diag=evidence['smd_spin_diagnostics'];matrix=np.load(source/'hessian_analytical_matrix.npy');logs=source/'analytical'
        work=p.with_suffix('');work.mkdir(parents=True)
        np.save(work/'hessian_cartesian.npy',matrix)
        for name in ['hessian','mol.xyz','thermo_25.log','thermo_50.log','thermo_100.log']:
            shutil.copy2(logs/name,work/name)
        spin=spin_degeneracy_eV(m)
        rec=dict(id=gid,state=state,status='derivative_checked_solution_thermal',identity=identity,converged=True,
            stationary_minimum=True,grid_response=True,smd_spin_diagnostics=spin_diag,
            gradient_max=evidence['gradient_max'],frequencies_cm=freq,derivative_checks=derivative,
            g_nuclear_qrrho_eV=values['50'],g_spin_partition_eV=spin,g_thermal_total_eV=values['50']+spin,
            nuclear_by_cutoff_eV=values,software=evidence['software'],implementation_controls=proof,
            adopted_from_verified_control=str(source),control_status_sha256=hashlib.sha256((source/'status.json').read_bytes()).hexdigest(),
            experimental_accuracy_validated=False,legacy_production_modified=False,created=time.time())
        temp=p.with_suffix('.tmp');temp.write_text(json.dumps(rec,indent=2));temp.replace(p)
        assert load_native_thermal(sd,q,m) is not None
        adoptions.append(dict(id=gid,state=state,native_record=str(p),control=str(source)))
    report=base/f'native_control_adoptions_{time.time_ns()}.json'
    report.write_text(json.dumps(dict(adopted=adoptions,legacy_production_modified=False),indent=2)+'\n')
    print(json.dumps(adoptions,indent=2))


if __name__=='__main__':main()
