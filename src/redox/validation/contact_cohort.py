"""Collect a declared chloride cohort without fixed outcome or duplicate counts."""
import argparse
import json
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.properties.contact_identity import validate_identities


def collect(roster_path, output):
    from redox.properties.conformer_geometry import alignment_molecule,symmetry_rmsd_matrix
    from redox.properties.conformer_ensemble import cluster_minima
    from redox.workflows.counterion_thermochemistry import checked_component
    if output.exists():raise FileExistsError('Preserve earlier contact cohorts')
    roster=json.loads(roster_path.read_text());starts=roster['declared_start_keys']
    for path,digest in roster['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Roster source changed')
    geometry_manifest=Path(roster['geometry_manifest'])
    if roster['source_sha256'].get(str(geometry_manifest))!=sha256_file(geometry_manifest):
        raise ValueError('Source-bound geometry roster required')
    if [r['key'] for r in json.loads(geometry_manifest.read_text())['states']]!=starts:
        raise ValueError('Declared starts differ from the original geometry campaign')
    identities=validate_identities(roster['records'],starts)
    rows=[];geometries=[];energies=[];context=None
    sources={str(roster_path):sha256_file(roster_path),str(Path(__file__)):sha256_file(__file__),**roster['source_sha256']}
    reviews=set()
    for item,(key,label) in zip(roster['records'],identities):
        review=Path(item['review']);ledger=Path(item['ledger'])
        if sha256_file(review)!=item['review_sha256'] or sha256_file(ledger)!=item['ledger_sha256']:
            raise ValueError('Declared component or cycle changed')
        if item['review_sha256'] in reviews:raise ValueError('Same checked component listed twice')
        reviews.add(item['review_sha256'])
        report,spec=checked_component(review)
        if report['key']!=item['original_start']:raise ValueError('Original start differs from reviewed identity')
        ledger_data=json.loads(ledger.read_text())
        if any(sha256_file(p)!=h for p,h in ledger_data['source_sha256'].items()):
            raise ValueError('Association ledger backing changed')
        cycle=next(r for r in ledger_data['states'] if r['key']==report['key'])
        if cycle['status']!='finite_primary_microstate_association_cycle_computed':
            raise ValueError('Completed matched association cycle required')
        for cutoff in ('25','50','100'):
            paired=cycle['component_contexts'][cutoff]['components']['paired']
            if (paired['sources']!={str(review):item['review_sha256']}
                    or paired['state']['smiles']!=report['smiles'] or paired['qc']!='independently_checked'):
                raise ValueError('Ledger must retain this exact checked contact component')
        current={k:report[k] for k in ('electronic_protocol','native_protocol','thermal_protocol',
                                      'temperature_K','standard_state','numerical_implementation')}
        if context is not None and current!=context:raise ValueError('Contact physical contexts differ')
        context=current;geometry=Path(spec['case'])/'state/opt.xyz'
        if sha256_file(geometry)!=report['geometry_sha256']:raise ValueError('Checked geometry changed')
        geometries.append(alignment_molecule(report['smiles'],geometry))
        energies.append(dict(key=key,E_eV=report['electronic_eV']))
        rows.append(dict(item,key=key,label=label,geometry=str(geometry),geometry_sha256=report['geometry_sha256'],
                         minimum_frequency_cm=report['minimum_frequency_cm'],
                         association_one_M_kJ_mol=cycle['one_M_deltaG_kJ_mol']))
        sources.update({str(p):sha256_file(p) for p in (review,ledger,geometry)})
    matrix,matches=symmetry_rmsd_matrix(geometries)
    groups={str(c):cluster_minima(energies,matrix,rmsd_cutoff_A=c,energy_tolerance_kJ_mol=.2) for c in (.05,.1,.2)}
    result=dict(status='declared_contact_cohort_duplicate_reviewed',identity_schema=2,records=rows,
        checked_components=len(rows),checked_outcomes=len(rows),declared_start_keys=starts,
        declared_original_starts=len(starts),distinct_original_starts_with_checked_G=len({r['original_start'] for r in rows}),
        context=context,proper_rotation_symmetry_heavy_RMSD_A=matrix.tolist(),automorphisms=matches,
        groups_by_RMSD_cutoff_A=groups,electronic_energy_tolerance_kJ_mol=.2,
        primary_association_range_kJ_mol=[min(r['association_one_M_kJ_mol']['50'] for r in rows),max(r['association_one_M_kJ_mol']['50'] for r in rows)],
        all_cutoff_associations_positive=all(v>0 for r in rows for v in r['association_one_M_kJ_mol'].values()),
        finite_contact_model_only=True,populations_inferred=False,conformer_ensemble_qualified=False,
        global_thermodynamic_pose_selected=False,material_admission=False,source_sha256=sources)
    for function in (validate_identities,alignment_molecule,symmetry_rmsd_matrix,cluster_minima):
        path=Path(function.__code__.co_filename);result['source_sha256'][str(path)]=sha256_file(path)
    output.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--roster',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=collect(args.roster.resolve(),args.out.resolve())
    print(json.dumps({k:result[k] for k in ('checked_outcomes','primary_association_range_kJ_mol','groups_by_RMSD_cutoff_A')}))
