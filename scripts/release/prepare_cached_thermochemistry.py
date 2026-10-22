#!/usr/bin/env python
"""Export a small, immutable raw-state arithmetic example, without launching QM."""
import argparse,hashlib,importlib.util,json,os,shutil
from pathlib import Path
from redox.core.common import ROOT,DFT,read_manifest
from redox.core.protocol import load_record,record_path
from redox.qm.native_thermal import cache_path,load_native_thermal
from redox.workflows.checkpoints import verify_state_checkpoint
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False)
    manifest={(r['id'],r['state']):r for r in read_manifest()}
    original=ROOT/'results/grafting/campaigns/2026-10-19_mecn_reference'
    base=ROOT/'results/grafting/campaigns/2026-10-19_base_assisted_pilot'
    states=[(name,DFT/gid/state,manifest[gid,state],None) for name,gid,state in [
        ('fc_neutral','ferrocene','neu'),('fc_oxidized','ferrocene','ox'),
        ('product_neutral','design_pmdi_ethyl','neu'),('product_reduced','design_pmdi_ethyl','red1')]]
    for name,key,campaign in [('precursor','gp_16aa455cb6d2',original),('site','merrifield_local_chloride',original),
            ('HCl','hydrogen_chloride',original),('chloride','chloride',original),('DBU','dbu_neutral',base),('DBUH','dbu_protonated',base)]:
        checked=verify_state_checkpoint(campaign/'done'/(key+'.json'))
        states.append((name,Path(checked['state_dir']),checked['task']['manifest_row'],checked['G_eV']))
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    exported=[];sources={}
    for name,sd,row,known in states:
        raw=json.loads((sd/'result.json').read_text());q=raw['charge'];mult=raw['mult'];energy=load_record(sd,q,mult)
        if energy is None:raise ValueError('Unqualified SP '+name)
        files={'geometry':'opt.xyz','result':'result.json','energy':str(record_path(sd,q,mult).relative_to(sd))}
        if (sd/'monatomic_thermal.json').exists():
            partition=json.loads((sd/'monatomic_thermal.json').read_text());thermal=partition['G_thermal_eV'];kind='monatomic'
            files['partition']='monatomic_thermal.json'
        else:
            partition=load_native_thermal(sd,q,mult)
            if partition is None:raise ValueError('Unqualified partition '+name)
            thermal=partition['g_thermal_total_eV'];kind='projected_DFT'
            native,_=cache_path(sd,q,mult)
            files.update(native=str(native.relative_to(sd)),partition=str(Path(partition['thermal_postprocessing_path']).relative_to(sd)),
                hessian=str((native.with_suffix('')/'hessian_cartesian.npy').relative_to(sd)))
        total=energy['e_smd_eV']+thermal
        if known is not None and abs(total-known)>1e-9:raise ValueError('Export differs from checked receipt '+name)
        dest=args.out/name;dest.mkdir();copied={}
        for role,rel in files.items():
            source=sd/rel;h=sha(source);sources[str(source.resolve())]=h
            target=dest/Path(rel).name
            if target.exists():raise ValueError('Export basename collision')
            shutil.copyfile(source,target)
            if sha(target)!=h:raise ValueError('Copy differs')
            copied[role]=str(target.relative_to(args.out))
        exported.append(dict(id=name,smiles=row['smiles'],charge=q,mult=mult,partition_kind=kind,files=copied,checked_G_eV=total))
    if any(sha(p)!=h for p,h in sources.items()):raise ValueError('A calculation changed during export')
    base_values=json.loads((ROOT/'results/grafting/base_assisted/current.json').read_text())
    expected=next(r for r in base_values['candidates'] if r['id']=='design_pmdi_ethyl')
    result=dict(schema=1,scope='Cached raw-state arithmetic, not a new QM calculation or external validation',
        producer_sha256=sha(__file__),states=exported,temperature_K=298.15,
        potential=dict(oxidized='product_neutral',reduced='product_reduced',reference_oxidized='fc_oxidized',reference_reduced='fc_neutral'),
        reactions=[dict(name='HCl_reference',coefficients={'product_neutral':1,'HCl':1,'precursor':-1,'site':-1},
            delta_n_solutes=0,source_kJ_mol=expected['reference_deltaG_kJ_mol']),
            dict(name='DBU_separated_ions',coefficients={'product_neutral':1,'DBUH':1,'chloride':1,'precursor':-1,'site':-1,'DBU':-1},
            delta_n_solutes=0,source_kJ_mol=expected['base_models']['dbu']['free_ions_kJ_mol'])],
        source_energy_conversion_kJ_mol_per_eV=96.485,SI_energy_conversion_kJ_mol_per_eV=96.48533212331002,
        original_sources_sha256=sources,files_sha256={str(p.relative_to(args.out)):sha(p) for p in args.out.rglob('*') if p.is_file()})
    (args.out/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(states=len(exported),files=len(result['files_sha256']),bytes=sum(p.stat().st_size for p in args.out.rglob('*') if p.is_file()))))


if __name__=='__main__':main()
