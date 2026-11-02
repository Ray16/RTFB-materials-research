"""Prepare frozen methyl-capped viologen motifs from checked polymer frames.

Reuses the archived constructor and its independent cap/pose auditor. Charged
viologens require an explicit CH2 wildcard: otherwise the symmetric core also
matches an N-methyl carbon, which is not the intended benzyl attachment atom.
This stage produces geometries, never interaction energies, G or populations.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.build.polymer_contacts import pi_core_groups, contact_metrics
from redox.workflows.molecular_native_namespace import verify, load

SOURCE_BYTES = Path(__file__).read_bytes()
KEYS = ('all_same__occ_2_3','all_same__occ_2_4','alternating__occ_2_3','alternating__occ_2_4')


def benzyl_query(fragment):
    if fragment.count('[*:1]') != 1:
        raise ValueError('One declared viologen attachment wildcard required')
    query = fragment.replace('[*:1]', '[*H2:1]')
    if Chem.MolFromSmarts(query) is None or Chem.MolFromSmiles(query) is None:
        raise ValueError('The constrained benzyl query is not supported')
    return query


def prepare(archive, preparation, stationary, output):
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT, archive)):
        raise FileExistsError('Fresh durable fragment namespace outside both source trees required')
    data = json.loads(preparation.read_text()); evaluated = json.loads(stationary.read_text())
    if (evaluated['status'] != 'exact_ionic_stationary_electronic_cycles_collected'
            or evaluated['native_G_qualified'] or evaluated['material_admission']):
        raise ValueError('The exact independently reviewed stationary ionic collection is required')
    cases = {r['case']['key']:r['case'] for r in data['cases']}
    query = benzyl_query(data['fragment'])
    graph_api = load(ROOT/'scripts/validation/review_charged_polymer_optimization.py', 'frozen_viologen_parent_mapping')
    construction_path = archive/'scripts/pipeline/local_motif_fragments.py'
    reviewer_path = archive/'scripts/validation/review_local_motif_fragments.py'
    constructor = load(construction_path, 'frozen_viologen_methyl_constructor')
    auditor = load(reviewer_path, 'frozen_viologen_methyl_auditor')
    rows, hashes = [], {str(p):sha256_file(p) for p in (preparation, stationary, construction_path, reviewer_path)}
    for key in KEYS:
        record = evaluated['records'][key]; path = Path(record['review'])
        if sha256_file(path) != record['review_sha256']:
            raise ValueError('Stationary component review changed')
        report = json.loads(path.read_text())
        if (report['status'] != 'charged_polymer_solution_stationary_geometry_independently_checked'
                or report['graph_key'] != key or (report['charge'],report['mult']) != (4,1)
                or report['polymer_smiles'] != cases[key]['smiles']
                or any(sha256_file(p) != h for p,h in report['source_sha256'].items())):
            raise ValueError('Exact charge-4/singlet parent frame and complete source review required')
        spec = verify(Path(report['namespace'])); geometry = Path(spec['sources'])/'input.xyz'
        if sha256_file(geometry) != report['geometry_sha256']:
            raise ValueError('Parent coordinate frame changed')
        parent = Chem.AddHs(Chem.MolFromSmiles(cases[key]['smiles']))
        observed = Chem.MolFromXYZFile(str(geometry)); rdDetermineBonds.DetermineConnectivity(observed,charge=0)
        mapping = graph_api._skeleton(observed).GetSubstructMatch(graph_api._skeleton(parent))
        if len(mapping) != parent.GetNumAtoms() or len(set(mapping)) != len(mapping):
            raise ValueError('A complete bijective graph-to-source coordinate map is required')
        xyz = read(geometry).positions[list(mapping)]
        # Independently check the returned permutation, including every hydrogen.
        expected_bonds = {tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in parent.GetBonds()}
        observed_bonds = {tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in observed.GetBonds()}
        mapped_bonds = {tuple(sorted((mapping[a],mapping[b]))) for a,b in expected_bonds}
        if (mapped_bonds != observed_bonds
                or [observed.GetAtomWithIdx(i).GetAtomicNum() for i in mapping] != [a.GetAtomicNum() for a in parent.GetAtoms()]):
            raise ValueError('The coordinate permutation does not preserve the complete parent graph')
        conformer = Chem.Conformer(len(mapping))
        for i,position in enumerate(xyz):conformer.SetAtomPosition(i,position)
        parent.AddConformer(conformer)
        groups = pi_core_groups(parent,query)
        parts = constructor.fragments(parent,query)
        if len(groups) != 2 or len(parts) != 2 or any((p['mol'].GetNumAtoms(),p['charge']) != (28,2) for p in parts):
            raise ValueError('Two exact methyl-viologen dications, each 28 atoms, are required')
        task = dict(cases[key], fragment=query, source_fragment=data['fragment'])
        rows.append(dict(key=key,task=task,atomic_numbers=[a.GetAtomicNum() for a in parent.GetAtoms()],
            coordinates_A=xyz.tolist(),groups=groups,contact=contact_metrics(xyz,groups)[0],
            source_geometry=str(geometry),source_geometry_sha256=sha256_file(geometry),
            graph_to_source_atom_indices=list(mapping)))
        hashes[str(path)] = sha256_file(path); hashes[str(geometry)] = sha256_file(geometry)
        hashes.update(report['source_sha256'])
    output.mkdir(parents=True)
    source = output/'sources'; source.mkdir()
    (source/'preparer.py').write_bytes(SOURCE_BYTES)
    hashes[str(source/'preparer.py')] = sha256_file(source/'preparer.py')
    feed = source/'spacing_conformers.json'
    atomic_write_json(feed,dict(states=rows,source_sha256=hashes,conformer_convergence_established=False,
        context=dict(solvent_model=evaluated['model'],quantum_settings=evaluated['quantum_settings'],
            numerical_implementation=evaluated['numerical_implementation']),
        formal_states='each capped motif +2 singlet; pair +4 singlet',
        scope='Four independently checked finite parent frames; no population or global minimum claim'))
    manifest = constructor.prepare(feed,output/'construction')
    checked = auditor.review(output/'construction',output/'construction_review.json')
    if checked['molecular_inputs'] != 36:
        raise ValueError('Complete four-frame, three-cap-length construction required')
    result = dict(status='four_ionic_parent_frames_methyl_fragment_construction_checked',
        construction=str(output/'construction/manifest.json'),construction_sha256=sha256_file(output/'construction/manifest.json'),
        independent_review=str(output/'construction_review.json'),independent_review_sha256=sha256_file(output/'construction_review.json'),
        source_frame_mapping=str(feed),source_frame_mapping_sha256=sha256_file(feed),
        capped_fragment_smiles=manifest['entries'][0]['molecules'][0]['smiles'],
        monomer_charge=2,pair_charge=4,multiplicities=1,frames=4,cap_lengths_A=[1.07,1.09,1.11],molecular_inputs=36,
        physical_context=json.loads(feed.read_text())['context'],
        interaction_definition='E_AB(frozen pair)-E_A(frozen)-E_B(frozen), each medium separately',
        source_cut='Benzyl CH2-to-polymer aryl bond replaced by C-H; retained benzyl carbon becomes methyl',
        cap_model='Dimethyl viologen dications; linker aryl/backbone removed, original motif coordinates retained',
        planned_energies='matched parent r2SCAN-D4/def2-SVP; gas and experimental smooth-surface SMD MeCN',
        interaction_energy_kJ_mol=None,free_energy_computed=False,quantum_jobs_submitted=0,
        deformation_energy_computed=False,basis_superposition_corrected=False,electronic_coupling_computed=False,
        conformer_population_inferred=False,full_polymer_nonadditivity_replaced=False,material_admission=False)
    atomic_write_json(output/'summary.json',result)
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','preparation','stationary','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare(a.archive.resolve(),a.preparation.resolve(),a.stationary.resolve(),a.out.resolve())
    print(json.dumps({k:result[k] for k in ('status','frames','molecular_inputs','monomer_charge','pair_charge')}))
