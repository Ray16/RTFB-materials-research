#!/usr/bin/env python
"""Compare preserved thermal models with matched Fc and conservative DFT reference.

An uncertainty envelope is model sensitivity, not a fitted Gaussian uncertainty.
One DFT state cannot validate thermal differences on all species or charge states.
"""
import argparse,csv,hashlib,json,math
import numpy as np
from pathlib import Path
from redox.core.common import ROOT,DFT,read_result
from redox.qm.thermal import paired_potential,spin_degeneracy_eV
from redox.qm.vibrations import parse_xtb_thermal_eV
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    if a.out.exists():raise ValueError('Choose a fresh analysis output; preserve previous evidence')
    a.out.mkdir(parents=True)
    records={};comparisons=[]
    for p in (BASE/'stationary').glob('*/*/*/result.json'):
        rec=json.loads(p.read_text());gid,st,model=rec['id'],rec['state'],rec['model']
        src=DFT/gid/st/'opt.xyz'
        current=src.exists() and hashlib.sha256(src.read_bytes()).hexdigest()==rec['source_geom_sha256']
        work=p.parent;values={}
        for k in (25,50,100):
            log=work/f'thermo_{k}.log'
            if log.exists():values[str(k)]=parse_xtb_thermal_eV(log.read_text())
        # SPH's mode-dependent bias scaling is not encoded by the saved Cartesian Hessian.
        # Raw-Hessian reprocessing cannot represent SPH at a different rotor cutoff.
        if model.endswith('sph'):values={'50':rec['g_thermal_eV']}
        valid=rec['n_imag']==0 and rec['min_frequency_cm']>=-5 and rec['connectivity_intact'] and rec['optimization_converged'] and model.endswith('relaxed')
        r=dict(id=gid,state=st,model=model,source_current=current,model_stationary_valid=valid,
            n_imag=rec['n_imag'],min_frequency_cm=rec['min_frequency_cm'],g_thermal_eV=rec['g_thermal_eV'],
            g_by_cutoff=values,source=str(p),transfer_validated=False)
        records[(gid,st,model)]=r;comparisons.append(r)
    # Independent mode-following can replace only this diagnostic model's failed minimum.
    for p in (BASE/'mode_follow_2026-10-11').glob('*/*/result.json'):
        rec=json.loads(p.read_text())
        if not rec['stationary_hessian_valid']:continue
        key=('aq_benzylamino',rec['state'],'alpb_minimum_mixed_accuracy_diagnostic')
        r=dict(id=key[0],state=key[1],model=key[2],source_current=True,model_stationary_valid=True,
            n_imag=rec['n_imag'],min_frequency_cm=rec['minimum_frequency_cm'],g_thermal_eV=rec['g_thermal_eV'],
            g_by_cutoff={},source=str(p),transfer_validated=False,total_energy_Ha=rec['total_energy_Ha'])
        if key not in records or r['total_energy_Ha']<records[key]['total_energy_Ha']:records[key]=r
    # For other states the same tight solution minimum definition uses stationary ALPB.
    for key,rec in list(records.items()):
        if rec['model']=='alpb_relaxed' and rec['model_stationary_valid']:
            k=(*key[:2],'alpb_minimum_mixed_accuracy_diagnostic')
            if k not in records:records[k]=dict(rec,model=k[2])
    matched_paths=list((BASE/'d3_correct_spin_solution_thermal_v2').glob('*/*/*/selected.json'))
    matched_paths+=list((BASE/'uniform_solution_thermal').glob('*/*/*/selected.json'))
    for p in matched_paths:
        rec=json.loads(p.read_text());sel=rec.get('selected')
        if not sel:continue
        r=dict(id=rec['id'],state=rec['state'],model='alpb_minimum_acc0p1',
            source_current=hashlib.sha256((DFT/rec['id']/rec['state']/'opt.xyz').read_bytes()).hexdigest()==rec['identity']['geometry_sha256'],
            model_stationary_valid=sel['stationary_valid'],n_imag=sel['n_imag'],min_frequency_cm=min(sel['frequencies_cm']),
            g_thermal_eV=sel['g_thermal_eV'],g_by_cutoff=rec['g_thermal_by_cutoff_eV'],source=str(p),transfer_validated=False,
            stable_branch_thermal_spread_eV=rec['stable_branch_thermal_spread_eV'])
        records[(r['id'],r['state'],r['model'])]=r;comparisons.append(r)
    model_rows=[]
    for gid in sorted({key[0] for key in records}):
        if gid=='ferrocene':continue
        names=sorted({key[1] for key in records if key[0]==gid},key=lambda st:-read_result(gid,st)['charge'])
        for model in sorted({key[2] for key in records if key[0]==gid}):
            keys=[(gid,st,model) for st in names]+[('ferrocene','ox',model),('ferrocene','neu',model)]
            allpresent=all(k in records and records[k]['source_current'] for k in keys)
            if not allpresent:continue
            stable=all(records[k]['model_stationary_valid'] for k in keys)
            for cutoff in (25,50,100):
                def g(gid,st):
                    r=records[(gid,st,model)]
                    return r['g_by_cutoff'].get(str(cutoff),r['g_thermal_eV'] if cutoff==50 else None)
                fcO,fcN=read_result('ferrocene','ox'),read_result('ferrocene','neu')
                for O,R in zip(names,names[1:]):
                    ro,rr=read_result(gid,O),read_result(gid,R)
                    if ro['charge']-rr['charge']!=1:continue
                    E=paired_potential(ro['e_smd_eV'],rr['e_smd_eV'],g(gid,O),g(gid,R),fcO['e_smd_eV'],fcN['e_smd_eV'],g('ferrocene','ox'),g('ferrocene','neu'))
                    if E is None:continue
                    model_rows.append(dict(id=gid,quantity=f'E({O}->{R})',model=model,cutoff_cm=cutoff,
                        value_V_or_eV=E,all_model_stationary=stable,
                        uniform_numerical_accuracy=model!='alpb_minimum_mixed_accuracy_diagnostic',transfer_validated=False))
                    spin_shift=(spin_degeneracy_eV(ro['mult'])-spin_degeneracy_eV(rr['mult'])
                        -spin_degeneracy_eV(fcO['mult'])+spin_degeneracy_eV(fcN['mult']))
                    model_rows[-1].update(spin_entropy_shift_V=spin_shift,
                        value_with_unsplit_spin_entropy=E+spin_shift,
                        electronic_entropy_assumption='unsplit spin multiplet only; no SOC/orbital degeneracy')
                if len(names)==3:
                    vals=[g(gid,s) for s in names]
                    if None not in vals:
                        G=[read_result(gid,s)['e_smd_eV']+v for s,v in zip(names,vals)]
                        model_rows.append(dict(id=gid,quantity=f'dG_disp({names[1]})',model=model,cutoff_cm=cutoff,
                            value_V_or_eV=G[0]+G[2]-2*G[1],all_model_stationary=stable,
                            uniform_numerical_accuracy=model!='alpb_minimum_mixed_accuracy_diagnostic',transfer_validated=False))
                        spin=[spin_degeneracy_eV(read_result(gid,s)['mult']) for s in names]
                        shift=spin[0]+spin[2]-2*spin[1]
                        model_rows[-1].update(spin_entropy_shift_V=shift,
                            value_with_unsplit_spin_entropy=G[0]+G[2]-2*G[1]+shift,
                            electronic_entropy_assumption='unsplit spin multiplet only; no SOC/orbital degeneracy')
    reference=BASE/'numerical_solution_hessian_gridresponse_verified/result.json'
    transfer=dict(reference_available=reference.exists(),production_modified=False,
        numerical_stationarity_is_not_transfer_accuracy=True,coverage_required='paired candidate charge states plus matched Fc/Fc+; calibrate transfer differences, not one absolute thermal correction')
    if reference.exists():
        ref=json.loads(reference.read_text())
        if ref['status']!='numerical_reference_verified' or not ref['stationary_minimum']:raise RuntimeError('Conservative DFT reference did not pass')
        transfer['reference']=str(reference);transfer['model_error_at_aq_red2_eV']={}
        for key,rec in records.items():
            if key[:2]!=('aq_benzylamino','red2') or not rec['source_current']:continue
            transfer['model_error_at_aq_red2_eV'][key[2]]=rec['g_thermal_eV']-ref['g_thermal_qrrho_by_cutoff_eV']['50']
        transfer['reference_g_thermal_eV']=ref['g_thermal_qrrho_by_cutoff_eV']['50']
        transfer['sorted_frequency_comparison']={}
        target=np.sort(np.asarray(ref['projected_frequencies_cm']))
        for key,rec in records.items():
            if key[:2]!=('aq_benzylamino','red2') or not rec['source_current']:continue
            source=json.loads(Path(rec['source']).read_text())
            frequencies=source.get('frequencies_cm')
            if frequencies is None and source.get('selected'):
                frequencies=source['selected']['frequencies_cm']
            if frequencies is None:continue
            low=np.asarray(frequencies)
            # Remove the six zero rigid modes before a sorted distribution comparison.
            if len(low)==len(target)+6:
                remove=np.argsort(abs(low))[:6];low=np.delete(low,remove)
            if len(low)!=len(target):continue
            low=np.sort(low);mask=target>100
            transfer['sorted_frequency_comparison'][key[2]]=dict(
                n_modes=len(low),high_frequency_mean_abs_relative_error=float(np.mean(abs(low[mask]-target[mask])/target[mask])),
                sorted_mean_absolute_error_cm=float(np.mean(abs(low-target))),
                note='Sorted distributions, not asserted one-to-one normal-mode identities.')
    def write(name,rows):
        if not rows:return
        flat=[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in rows]
        fields=list(dict.fromkeys(k for r in flat for k in r))
        with (a.out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(flat)
    write('state_thermal_models.csv',comparisons);write('matched_fc_predictions.csv',model_rows)
    (a.out/'transfer_reference.json').write_text(json.dumps(transfer,indent=2))
    print(json.dumps(dict(states_compared=len(comparisons),paired_predictions=len(model_rows),reference=transfer),indent=2))


if __name__=='__main__':main()
