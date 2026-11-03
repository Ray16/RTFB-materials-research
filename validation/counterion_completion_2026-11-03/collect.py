"""Collect six checked outcomes without counting repeated searches as populations."""
import json
import os
from pathlib import Path

from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.properties.conformer_geometry import alignment_molecule,symmetry_rmsd_matrix
from redox.properties.conformer_ensemble import cluster_minima
from redox.workflows.counterion_thermochemistry import checked_component

OUT=Path(__file__).resolve().parent
previous=ROOT/'validation/counterion_sensitivity_2026-11-02/thermochemistry_current_04.json'
new_ledger=OUT/'association_new_contacts.json'
old=json.loads(previous.read_text());new=json.loads(new_ledger.read_text())
rows=[];sources={str(p):sha256_file(p) for p in (previous,new_ledger,Path(__file__))}
items=[dict(r) for r in old['records']]
for key in ('N0_below','N0_outer'):
    value=next(r for r in new['states'] if r['key']==key)
    items.append(dict(original_start=key,branch='minus',review=str(OUT/(key+'_minus_review.json')),
        ledger=str(new_ledger),association_one_M_kJ_mol=value['one_M_deltaG_kJ_mol']))
geometries=[];energy_rows=[];common=None
for item in items:
    path=Path(item['review']);r,spec=checked_component(path)
    if item['original_start']!=r['key']:raise ValueError('Outcome label differs from its original start')
    ledger=Path(item['ledger']);cycle=next(c for c in json.loads(ledger.read_text())['states'] if c['key']==r['key'])
    if item['association_one_M_kJ_mol']!=cycle['one_M_deltaG_kJ_mol']:raise ValueError('Association values differ')
    context={k:r[k] for k in ('electronic_protocol','native_protocol','thermal_protocol','temperature_K','standard_state','numerical_implementation')}
    if common is not None and context!=common:raise ValueError('Unmatched component contexts')
    common=context;geometry=Path(spec['case'])/'state/opt.xyz'
    if sha256_file(geometry)!=r['geometry_sha256']:raise ValueError('Geometry differs from checked thermochemistry')
    label=r['key']+('_'+item['branch'] if item['branch'] else '')
    geometries.append(alignment_molecule(r['smiles'],geometry));energy_rows.append(dict(key=label,E_eV=r['electronic_eV']))
    rows.append(dict(key=label,original_start=r['key'],branch=item['branch'],review=str(path),
        review_sha256=sha256_file(path),ledger=str(ledger),ledger_sha256=sha256_file(ledger),
        minimum_frequency_cm=r['minimum_frequency_cm'],association_one_M_kJ_mol=cycle['one_M_deltaG_kJ_mol']))
    sources.update({str(p):sha256_file(p) for p in (path,ledger,geometry)})
matrix,matches=symmetry_rmsd_matrix(geometries)
groups={str(c):cluster_minima(energy_rows,matrix,rmsd_cutoff_A=c,energy_tolerance_kJ_mol=.2) for c in (.05,.1,.2)}
result=dict(status='six_checked_contact_outcomes_duplicate_reviewed',records=rows,
    checked_outcomes=len(rows),distinct_original_starts_with_checked_G=len({r['original_start'] for r in rows}),
    context=common,proper_rotation_symmetry_heavy_RMSD_A=matrix.tolist(),automorphisms=matches,
    groups_by_RMSD_cutoff_A=groups,electronic_energy_tolerance_kJ_mol=.2,
    primary_association_range_kJ_mol=[min(r['association_one_M_kJ_mol']['50'] for r in rows),max(r['association_one_M_kJ_mol']['50'] for r in rows)],
    all_primary_associations_positive=all(r['association_one_M_kJ_mol']['50']>0 for r in rows),
    finite_contact_model_only=True,populations_inferred=False,conformer_ensemble_qualified=False,
    material_admission=False,manuscript_frozen_values_replaced=False,source_sha256=sources)
for name in ('conformer_geometry.py','conformer_ensemble.py'):
    p=ROOT/'src/redox/properties'/name;result['source_sha256'][str(p)]=sha256_file(p)
path=OUT/'current.json'
if path.exists():raise FileExistsError('Preserve prior collections')
path.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('checked_outcomes','primary_association_range_kJ_mol','groups_by_RMSD_cutoff_A')}))
