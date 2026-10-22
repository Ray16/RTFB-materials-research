#!/usr/bin/env python3
"""Recompute one potential and grafting sums from portable cached state records.

Only the Python standard library is used. Original absolute paths are descriptive
metadata; no original source location, quantum package or xTB library is opened.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


def reconstruct(folder):
    folder=Path(folder).resolve();manifest=json.loads((folder/'manifest.json').read_text())
    def path(relative):
        p=(folder/relative).resolve()
        if not p.is_relative_to(folder):raise ValueError('Cached path escapes the example')
        return p
    for relative,expected in manifest['files_sha256'].items():
        if hashlib.sha256(path(relative).read_bytes()).hexdigest()!=expected:
            raise ValueError('Cached record hash mismatch: '+relative)
    def load(relative):return json.loads(path(relative).read_text())
    states={};inventories={};charges={}
    for row in manifest['states']:
        files=row['files'];geo=path(files['geometry']).read_bytes();lines=geo.decode().splitlines();n=int(lines[0])
        atoms=[v.split() for v in lines[2:]]
        if len(atoms)!=n or any(len(a)!=4 or not all(math.isfinite(float(v)) for v in a[1:]) for a in atoms):
            raise ValueError('Malformed cached geometry')
        raw=load(files['result']);energy=load(files['energy']);partition=load(files['partition'])
        if any(d['charge']!=row['charge'] or d['mult']!=row['mult'] for d in (raw,energy)):
            raise ValueError('Charge/spin mismatch')
        if energy['geom_sha1']!=hashlib.sha1(geo).hexdigest()[:16] or not energy['converged_smd']:
            raise ValueError('Electronic energy does not match the cached geometry')
        if row['partition_kind']=='monatomic':
            if n!=1 or partition['internal_modes']!=0:raise ValueError('Not monatomic')
            thermal=partition['G_thermal_eV']
        else:
            native=load(files['native']);gh=hashlib.sha256(geo).hexdigest()
            for d in (native,partition):
                if d['identity']['geometry_sha256']!=gh or d['identity']['charge']!=row['charge'] or d['identity']['mult']!=row['mult']:
                    raise ValueError('Thermal identity mismatch')
            if partition['identity']['raw_native_sha256']!=hashlib.sha256(path(files['native']).read_bytes()).hexdigest():
                raise ValueError('Partition/native mismatch')
            modes=partition['frequencies_cm']
            if len(modes)!=partition['internal_mode_count'] or len(modes)!=3*n-partition['rigid_rank'] or min(modes)<-5:
                raise ValueError('Invalid cached mode list')
            spin=-8.617333262145e-5*manifest['temperature_K']*math.log(row['mult'])
            if abs(spin-partition['g_spin_partition_eV'])>1e-8:raise ValueError('Spin partition differs')
            thermal=partition['nuclear_by_cutoff_eV']['50']+spin
            if abs(thermal-partition['g_thermal_total_eV'])>1e-8:raise ValueError('Cached partition sum differs')
        value=energy['e_smd_eV']+thermal
        if not math.isfinite(value) or abs(value-row['checked_G_eV'])>1e-8:raise ValueError('State G differs')
        states[row['id']]=value;inventories[row['id']]=Counter(a[0] for a in atoms);charges[row['id']]=row['charge']
    p=manifest['potential'];fc=states[p['reference_oxidized']]-states[p['reference_reduced']]
    potential=states[p['oxidized']]-states[p['reduced']]-fc
    reactions=[]
    for reaction in manifest['reactions']:
        coefficients=reaction['coefficients'];balance=Counter()
        for key,nu in coefficients.items():
            for atom,count in inventories[key].items():balance[atom]+=nu*count
        if any(balance.values()) or sum(charges[key]*nu for key,nu in coefficients.items()):
            raise ValueError('Reaction atom/charge balance fails')
        if sum(coefficients.values())!=reaction['delta_n_solutes']:raise ValueError('Solute-count mismatch')
        dg=math.fsum(states[key]*nu for key,nu in coefficients.items())
        standard=reaction['delta_n_solutes']*.00831446261815324*manifest['temperature_K']*math.log(.08205736608096*manifest['temperature_K'])
        legacy=dg*manifest['source_energy_conversion_kJ_mol_per_eV']+standard
        if abs(legacy-reaction['source_kJ_mol'])>1e-6:raise ValueError('Reaction differs from checked source')
        reactions.append(dict(name=reaction['name'],deltaG_eV=dg,SI_deltaG_kJ_mol=dg*manifest['SI_energy_conversion_kJ_mol_per_eV']+standard,
            source_factor_deltaG_kJ_mol=legacy,standard_state_correction_kJ_mol=standard))
    return dict(status='cached_state_arithmetic_reproduced',states=len(states),Fc_absolute_V=fc,
        PMDI_first_reduction_V_vs_Fc=potential,reactions=reactions,
        quantum_calculations_run=False,nuclear_partition_reintegrated=False,external_accuracy_inferred=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path(__file__).resolve().parent/'data')
    parser.add_argument('--out',type=Path);args=parser.parse_args();result=reconstruct(args.data)
    if args.out:
        with args.out.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
