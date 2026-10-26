#!/usr/bin/env python
"""Independently reconstruct the exact local TEMPO oxidation/Fc comparison."""
import argparse
import json
import math
from pathlib import Path

from rdkit import Chem

from redox.core.common import ROOT, read_result
from redox.core.protocol import load_record, record_path
from redox.qm.native_thermal import cache_path
from redox.workflows import polymer_campaign as api


def potential(oxidized, reduced, fc_oxidized, fc_reduced):
    values = (oxidized,reduced,fc_oxidized,fc_reduced)
    if any(type(v) not in (int,float) or not math.isfinite(v) for v in values):
        raise ValueError('Finite native state free energies required')
    return math.fsum((oxidized,-reduced,-fc_oxidized,fc_reduced))


def review(route_path, out):
    route = json.loads(route_path.read_text())
    product = route['candidate']['product_smiles']
    if route['candidate']['id'] != 'tempo' or Chem.MolToInchiKey(Chem.MolFromSmiles(product)) == '':
        raise ValueError('Expected the exact TEMPO local graft')
    with (ROOT/'library/manifest.csv').open() as stream:
        import csv
        matches = [r for r in csv.DictReader(stream) if r['id'] == 'tempo' and r['state'] in ('rad','ox')]
    if len(matches) != 2 or any(Chem.MolToSmiles(Chem.MolFromSmiles(r['smiles'])) != product for r in matches):
        raise ValueError('TEMPO graph differs from the active state manifest')
    audit = api.load_module('independent_tempo_native_partition',ROOT/'scripts/validation/audit/native_partition.py')
    sources = {str(p):api.digest(p) for p in (route_path,Path(__file__).resolve(),
               ROOT/'scripts/validation/audit/native_partition.py',ROOT/'library/manifest.csv')}
    states = {}
    expected = [('tempo','rad',0,2),('tempo','ox',1,1),('ferrocene','neu',0,1),('ferrocene','ox',1,2)]
    for gid,state,q,mult in expected:
        directory = ROOT/'calcs/dft'/gid/state
        raw,selected = read_result(gid,state,raw=True),read_result(gid,state)
        if raw is None or selected is None or (raw['charge'],raw['mult']) != (q,mult) or selected.get('geometry_spin_qc','ok') != 'ok':
            raise ValueError('Geometry, state or spin gate failed')
        energy = load_record(directory,q,mult)
        surface,_ = cache_path(directory,q,mult)
        partition = audit.read_partition(surface,directory/'opt.xyz')
        if energy is None or partition is None or energy.get('smd_spin_contam',1) > 0.1:
            raise ValueError('Active energy or independent native partition failed')
        if selected.get('thermal_qc') != 'ok' or selected.get('spin_ensemble_status') not in ('not_selected',None):
            raise ValueError('Native thermal selection is not admitted')
        total = energy['e_smd_eV']+partition['g_thermal_total_eV']
        by_cutoff = {str(c):energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][str(c)]
                     +partition['g_thermal_total_eV']-partition['nuclear_by_cutoff_eV']['50'] for c in (25,50,100)}
        if abs(total-by_cutoff['50']) > 1e-8:
            raise ValueError('Independent thermal reconstruction differs')
        record_file=record_path(directory,q,mult)
        for path in (directory/'result.json',directory/'opt.xyz',record_file,surface):
            sources[str(path)] = api.digest(path)
        sources.update(partition['evidence_sha256'])
        states[gid+'__'+state] = dict(charge=q,multiplicity=mult,active_E_smd_eV=energy['e_smd_eV'],
            native_thermal_eV=partition['g_thermal_total_eV'],G_eV=total,G_by_cutoff_eV=by_cutoff,
            smd_spin_contamination=energy['smd_spin_contam'],partition_path=partition['path'])
    values = {str(c):potential(states['tempo__ox']['G_by_cutoff_eV'][str(c)],
              states['tempo__rad']['G_by_cutoff_eV'][str(c)],states['ferrocene__ox']['G_by_cutoff_eV'][str(c)],
              states['ferrocene__neu']['G_by_cutoff_eV'][str(c)]) for c in (25,50,100)}
    api.verify_sources(sources)
    result=dict(status='TEMPO_local_oxidation_native_partition_independently_reconstructed',
        product_smiles=product,product_inchikey=Chem.MolToInchiKey(Chem.MolFromSmiles(product)),
        event='ox->rad',potential_V_vs_Fc_by_rotor_cutoff=values,states=states,source_sha256=sources,
        experimental_exact_graft_potential_claimed=False,free_energy_of_grafting_computed=False,
        oligomer_interaction_or_loading_inferred=False,all_axis_eligibility=False,
        scope='Unfitted local graft-site MeCN model; one-electron nitroxide/oxoammonium comparison; finite conformer and Fc-reference limitations remain')
    if out.exists():raise FileExistsError('Preserve an existing review')
    out.write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--route',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=review(a.route.resolve(),a.out.resolve());print(json.dumps({'status':r['status'],'potential_V':r['potential_V_vs_Fc_by_rotor_cutoff']['50']}))
