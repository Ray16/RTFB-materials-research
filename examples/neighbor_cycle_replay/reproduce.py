#!/usr/bin/env python3
"""Replay the exact-PMDI finite neighbor cycle from contained cached records."""
import argparse
import importlib.util
import json
import math
from pathlib import Path


def native_replay(root,manifest):
    """Independently reconstruct modes and nuclear sums; no quantum library."""
    import numpy as np
    modules={}
    for name in ('projected_modes','nuclear_equations'):
        spec=importlib.util.spec_from_file_location('neighbor_'+name,root/(name+'.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);modules[name]=module
    checked=[];monatomic=[]
    def local(relative):
        path=(root/relative).resolve()
        if not path.is_relative_to(root):raise ValueError('Native artifact escapes package')
        return path
    for row in manifest['states']:
        files=row['files'];part=json.loads(local(files['partition']).read_text())
        lines=local(files['geometry']).read_text().splitlines()[2:]
        xyz=np.array([[float(v) for v in line.split()[1:]] for line in lines])
        if row['partition_kind']=='monatomic':
            mass=float(part['mass_amu'])
            if not math.isfinite(mass) or mass<=0:raise ValueError('Invalid atomic mass')
            value=modules['nuclear_equations'].nuclear_equations(np.array([]),xyz,np.array([mass]),1,50)
            error=abs(value-part['G_thermal_eV'])
            if not math.isfinite(error) or error>1e-8:raise ValueError('Monatomic partition differs')
            monatomic.append(dict(id=row['id'],nuclear_error_eV=error));continue
        mass=np.asarray(part['atomic_masses_amu'])
        hessian=np.load(local(files['hessian']),allow_pickle=False)
        if (mass.shape!=(len(xyz),) or not np.isfinite(mass).all() or np.any(mass<=0)
                or hessian.shape!=(3*len(xyz),3*len(xyz)) or not np.isfinite(hessian).all()
                or abs(hessian-hessian.T).max()>1e-6):raise ValueError('Invalid cached Hessian or masses')
        frequencies,rank=modules['projected_modes'].projected_modes(hessian,xyz,mass)
        expected=np.asarray(part['frequencies_cm'])
        if (expected.shape!=frequencies.shape or not np.isfinite(expected).all()
                or len(frequencies)!=part['internal_mode_count'] or rank!=part['rigid_rank']
                or min(frequencies)<-5):raise ValueError('Invalid native internal-mode coverage/minimum')
        error=float(np.max(abs(frequencies-expected)))
        if error>1e-5:raise ValueError('Independent Hessian modes differ')
        sigma=part['symmetry_number']
        if type(sigma) is not int or sigma<1:raise ValueError('Positive integral symmetry number required')
        thermal={}
        for cutoff in (25,50,100):
            value=modules['nuclear_equations'].nuclear_equations(frequencies,xyz,mass,sigma,cutoff)
            difference=abs(value-part['nuclear_by_cutoff_eV'][str(cutoff)])
            if not math.isfinite(difference) or difference>1e-8:raise ValueError('Independent native nuclear sum differs')
            thermal[str(cutoff)]=difference
        checked.append(dict(id=row['id'],internal_modes=len(frequencies),symmetry_number=sigma,
            frequency_max_error_cm=error,nuclear_errors_eV=thermal))
    return dict(status='raw_neighbor_Hessians_modes_and_nuclear_sums_independently_reconstructed',
        Hessians_reconstructed=len(checked),monatomic_partitions_recomputed=len(monatomic),
        states=checked,monatomic_states=monatomic,numpy_version=np.__version__,
        original_stationarity_spin_and_displaced_SCF_QC_repeated=False,xTB_library_loaded=False)


def replay(root,full_native=False):
    root=Path(root).resolve()
    manifest=json.loads((root/'manifest.json').read_text())
    if (manifest['synthesis_admission'] or manifest['conformer_convergence_established']
            or manifest['isolated_motif_interaction_computed'] or manifest['electronic_coupling_computed']):
        raise ValueError('Cached finite-cycle replay cannot promote its scientific scope')
    if any(type(row['checked_G_eV']) not in (int,float) or not math.isfinite(row['checked_G_eV'])
           for row in manifest['states']):raise ValueError('Finite checked state G required')
    names=[row['id'] for row in manifest['states']]
    if len(names)!=14 or len(set(names))!=14:raise ValueError('Fourteen distinct raw state identities required')
    def number(value):
        if type(value) not in (int,float) or not math.isfinite(value):raise ValueError('Finite source value required')
        return value
    spec=importlib.util.spec_from_file_location('contained_cached_state_reader',root/'cached_reader.py')
    reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
    primary=reader.reconstruct(root)
    def load(relative):
        path=(root/relative).resolve()
        if not path.is_relative_to(root):raise ValueError('Cached artifact escapes package')
        return json.loads(path.read_text())
    source=load(manifest['neighbor_source'])
    if (source['candidate_id']!='design_pmdi_ethyl' or source['solvent']!='acetonitrile'
            or source['conformer_convergence_established'] or source['synthesis_admission']):
        raise ValueError('Expected finite exact-PMDI cycle source required')
    aliases={'pair':'polymer_pair','empty':'empty','site2':'site2','site3':'site3'}
    states={row['id']:row for row in manifest['states']}
    conversion=number(manifest['SI_energy_conversion_kJ_mol_per_eV'])
    if conversion<=0 or number(manifest['temperature_K'])!=298.15:raise ValueError('Expected positive conversion and 298.15 K')
    standard=-.00831446261815324*manifest['temperature_K']*math.log(.08205736608096*manifest['temperature_K'])
    cutoffs={}
    for cutoff in ('25','50','100'):
        values={}
        for key,row in states.items():
            energy=load(row['files']['energy']);nuclear=load(row['files']['partition'])
            correction=(nuclear['G_thermal_eV'] if row['partition_kind']=='monatomic'
                else nuclear['nuclear_by_cutoff_eV'][cutoff]+nuclear['g_spin_partition_eV'])
            value=energy['e_smd_eV']+correction
            if not math.isfinite(value):raise ValueError('Nonfinite cutoff-specific G')
            values[key]=value
        expected=source['by_cutoff'][cutoff]
        for name,key in aliases.items():
            if abs(values[key]-number(expected['state_G_eV'][name]))>1e-8:raise ValueError('Cycle state differs')
        delta=math.fsum([values['polymer_pair'],values['empty'],-values['site2'],-values['site3']])*conversion
        if abs(delta-number(expected['deltaDeltaG_kJ_mol']))>1e-6:raise ValueError('Neighbor nonadditivity differs')
        models={}
        for mode in ('free_ions','contact_pair'):
            reagent=math.fsum([values['protonated'],values['chloride'],-values['precursor'],-values['base']]) if mode=='free_ions' else math.fsum([values['pair'],-values['precursor'],-values['base']])
            correction=0. if mode=='free_ions' else standard
            computed={}
            for site,other in ((2,3),(3,2)):
                computed[f'graft_{site}_on_empty_kJ_mol']=math.fsum([values[f'site{site}'],-values['empty'],reagent])*conversion+correction
                computed[f'graft_{site}_given_{other}_kJ_mol']=math.fsum([values['polymer_pair'],-values[f'site{other}'],reagent])*conversion+correction
            computed['direct_double_graft_kJ_mol']=math.fsum([values['polymer_pair'],-values['empty'],2*reagent])*conversion+2*correction
            for name,value in computed.items():
                if abs(value-number(expected['models'][mode][name]))>1e-6:raise ValueError('Conditional reaction differs')
            for site,other in ((2,3),(3,2)):
                closure=computed[f'graft_{site}_on_empty_kJ_mol']+computed[f'graft_{other}_given_{site}_kJ_mol']
                if abs(closure-computed['direct_double_graft_kJ_mol'])>1e-6:raise ValueError('Attachment-order cycle does not close')
            models[mode]=computed
        cutoffs[cutoff]=dict(deltaDeltaG_kJ_mol=delta,models=models)
    result=dict(status='portable_exact_PMDI_neighbor_cycle_reconstructed',states=primary['states'],
        balanced_reactions=len(primary['reactions']),PMDI_first_reduction_V_vs_Fc=primary['PMDI_first_reduction_V_vs_Fc'],
        by_cutoff=cutoffs,quantum_calculations_run=False,nuclear_partition_reintegrated=False,
        conformer_convergence_established=False,synthesis_admission=False,
        isolated_motif_interaction_computed=False,electronic_coupling_computed=False)
    if full_native:
        result['native_reconstruction']=native_replay(root,manifest)
        result['nuclear_partition_reintegrated']=True
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument('--out',type=Path);parser.add_argument('--full-native',action='store_true')
    args=parser.parse_args();result=replay(args.data,args.full_native)
    if args.out:
        with args.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
