#!/usr/bin/env python3
"""Compare two scoped electronic quantities without relabeling their difference."""
import argparse
import json
import math
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.motif_interaction_energy import checked_triplet
from redox.validation.execution_review import freeze_source

SOURCE_BYTES=Path(__file__).read_bytes()


def compare(collection,stationary,output):
    if output.exists():raise FileExistsError('Preserve earlier interaction-definition comparisons')
    local=json.loads(collection.read_text());whole=json.loads(stationary.read_text())
    if local['status']!='twelve_triplet_frozen_charged_motif_E_collection' or local['independent_review_errors']:
        raise ValueError('A complete-roster reviewed local-motif ledger without errors is required')
    if whole['status']!='exact_ionic_stationary_electronic_cycles_collected':raise ValueError('Exact balanced whole-chain E cycles required')
    for path,digest in whole['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Whole-chain electronic evidence changed')
    hashes={str(p):sha256_file(p) for p in (collection,stationary,freeze_source(output.parent,Path(__file__),SOURCE_BYTES))}
    values={};reported=0
    for row in local['states']:
        if row['interaction_E_kJ_mol_by_medium'] is None:continue
        path=Path(row['review'])
        if sha256_file(path)!=row['review_sha256']:raise ValueError('Local electronic review changed')
        checked,_=checked_triplet(path,row['key'])
        if checked['interaction_E_kJ_mol_by_medium']!=row['interaction_E_kJ_mol_by_medium']:
            raise ValueError('Local interaction ledger differs from its raw-backed review')
        expected=whole['quantum_settings'];method=checked['method']
        keys={'xc':'xc','basis':'basis','disp':'disp','grid_level':'grid_level',
              'conv_tol':'SCF_conv_tol','conv_tol_grad':'SCF_conv_tol_grad','max_cycle':'SCF_max_cycle'}
        if (checked['solvent_model']!=whole['model'] or checked['numerical_implementation']!=whole['numerical_implementation']
                or any(method[key]!=expected[parent_key] for key,parent_key in keys.items())):
            raise ValueError('Physical/numerical models must match before comparing the definitions')
        hashes[str(path)]=sha256_file(path)
        if checked['cap_length_A']==1.09:values[checked['parent_key']]=checked
    rows=[]
    for cycle in whole['cycles']:
        definition=cycle['definition'];key=definition['configuration']+'__occ_'+'_'.join(map(str,definition['sites']))
        matched=values.get(key);full=cycle['delta_delta_E_kJ_mol']
        if full is None or not math.isfinite(full):raise ValueError('The declared whole-chain electronic cycle must be complete')
        motif=None if matched is None else matched['interaction_E_kJ_mol_by_medium']['acetonitrile']
        rows.append(dict(parent_key=key,whole_chain_nonadditivity_E_kJ_mol=full,
            frozen_methyl_motif_interaction_E_kJ_mol=motif,
            same_sign=None if motif is None else (full>=0)==(motif>=0),
            cap_sensitivity=local['cap_sensitivity'].get(key),
            unique_backbone_or_deformation_energy_assigned=False))
    result=dict(status='scoped_whole_chain_and_frozen_motif_E_definitions_compared',states=rows,
        compared_central_cap_frames=sum(r['frozen_methyl_motif_interaction_E_kJ_mol'] is not None for r in rows),
        opposite_sign_frames=sum(r['same_sign'] is False for r in rows),
        definition_whole='E(Pij)+E(P0)-E(Pi)-E(Pj), with each parent graph at its retained optimized geometry',
        definition_local='E(MA+MB)-E(MA)-E(MB), with methyl caps and parent-relative frozen motif coordinates',
        interpretation='A favourable whole-chain occupancy cycle is not an isolated motif-attraction measurement. The quantities differ in backbone/linker/cap, deformation and continuum contributions; their difference is not assigned uniquely.',
        free_energies_compared=False,electronic_coupling_computed=False,experimental_validation_inferred=False,
        material_admission=False,source_sha256=hashes)
    atomic_write_json(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('collection','stationary','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=compare(a.collection.resolve(),a.stationary.resolve(),a.out.resolve())
    print(json.dumps(dict(status=r['status'],compared_frames=r['compared_central_cap_frames'],opposite_sign_frames=r['opposite_sign_frames'])))
