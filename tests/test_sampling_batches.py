"""Finite cells retain source/context checks without inventing missing cycles."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from redox.core.hashing import sha256_file
from redox.validation import sampling_batches as api


def fixture(tmp_path, monkeypatch):
    contract = json.loads((Path(__file__).resolve().parents[1]/
        'validation/corrected_sampling_2026-11-02/sampling_contract.json').read_text())
    graph = {k: [] for k in api.GRAPH_CONTEXT}
    graph.update(smiles='CC', mapped_smiles='CC', graph_sha256='graph', model_key='all_same__occ_2_3',
                 charge=0, mult=1, atoms=8, chain_ends=2, backbone_carbon_count=2, backbone_CC_bonds=1)
    context = {k: contract[k] for k in ('experimental_method_id', 'electronic_protocol', 'native_settings')}
    context.update(experimental_model=contract['solvent_model'], charge=0, multiplicity=1,
        temperature_K=298.15, standard_state='declared translation',
        thermal_protocol=dict(pressure='declared translation', temperature_K=298.15, rotor_cutoffs_cm=[25, 50, 100]))
    rows = []; records = []; sources = {}
    for configuration in contract['configurations']:
        for sites in contract['occupancy_graphs']:
            label = configuration+'__occ_'+('_'.join(map(str, sites)) if sites else 'none')
            for i in contract['expanded_batch_input_IDs']:
                key = label+f'__conf{i:03d}'
                row = dict(key=key, status='pending', native_G_qualified=False,
                           G_by_cutoff_eV=None, nuclear_by_cutoff_eV=None)
                if label == 'all_same__occ_2_3':
                    folder = tmp_path/key; folder.mkdir()
                    geometry = folder/'input.xyz'; geometry.write_text('fixture coordinates')
                    namespace = folder/'namespace.json'; namespace.write_text(json.dumps(dict(case=str(folder))))
                    review = folder/'review.json'; review.write_text('checked fixture receipt')
                    sources[str(review)] = sha256_file(review)
                    ctx = deepcopy(context); ctx.update(task=dict(graph, key=key), geometry_sha256=sha256_file(geometry))
                    row.update(status='native_independently_checked', native_G_qualified=True, context=ctx,
                        namespace=str(namespace), review=str(review), review_sha256=sources[str(review)],
                        electronic_E_smd_eV=-100., nuclear_by_cutoff_eV={str(k): 1. for k in (25, 50, 100)},
                        G_by_cutoff_eV={str(k): -99. for k in (25, 50, 100)})
                    records.append(dict(key=key, review=str(review), G_by_cutoff_eV=row['G_by_cutoff_eV']))
                rows.append(row)
    native = dict(status='complete_roster_corrected_native_review_snapshot', states=rows,
        declared_inputs=240, independent_review_errors=[], occupied_ensemble_cycles_qualified=False,
        pending_values_are_null=True, checked_native_inputs=20)
    audit = dict(status='frozen_sampling_component_receipts_rechecked', records=records,
        producer_sha256=sha256_file(Path(api.__file__).with_name('sampling_snapshot.py')),
        context_source_and_component_sums_rechecked=True, occupied_ensemble_cycles_qualified=False,
        primary_candidate_admission=False, checked_native_inputs=20, source_sha256=sources)
    paths = [tmp_path/name for name in ('native.json', 'audit.json', 'contract.json')]
    monkeypatch.setattr(api, 'alignment_molecule', lambda *a: object())
    monkeypatch.setattr(api, 'symmetry_rmsd_matrix', lambda mols: (np.ones((len(mols), len(mols)))-np.eye(len(mols)), 1))
    return paths, native, audit, contract


def run(paths, native, audit, contract):
    paths[0].write_text(json.dumps(native))
    audit['source_sha256'][str(paths[0])] = sha256_file(paths[0])
    paths[1].write_text(json.dumps(audit)); paths[2].write_text(json.dumps(contract))
    return api.analyze(*paths)


def test_complete_cell_is_compared_but_missing_cells_and_cycles_stay_unqualified(tmp_path, monkeypatch):
    args = fixture(tmp_path, monkeypatch); result = run(*args)
    assert result['complete_cells'] == 1
    assert len(result['cells']) == 12
    assert sum(c['partitions'] is None for c in result['cells']) == 11
    assert result['occupied_ensemble_cycles_qualified'] is False
    cell = next(c for c in result['cells'] if c['partitions'] is not None)
    assert cell['partitions']['0.1']['50']['expanded']['unique_minima'] == 20
    assert cell['partitions']['0.1']['50']['comparison']['ensemble_G_change_kJ_mol'] == pytest.approx(-1.7182820758)


@pytest.mark.parametrize('bad', ['solvent', 'charge', 'spin', 'temperature', 'pressure', 'thermal',
                                 'graph', 'component', 'geometry', 'receipt', 'roster', 'audit', 'batches'])
def test_mixed_context_changed_sources_and_invented_inputs_are_refused(tmp_path, monkeypatch, bad):
    paths, native, audit, contract = fixture(tmp_path, monkeypatch)
    row = next(r for r in native['states'] if r['native_G_qualified']); ctx = row['context']
    if bad == 'solvent': ctx['experimental_model']['solvent'] = 'water'
    if bad == 'charge': ctx['charge'] = 1
    if bad == 'spin': ctx['multiplicity'] = 3
    if bad == 'temperature': ctx['temperature_K'] = 300
    if bad == 'pressure': ctx['standard_state'] = 'different pressure'
    if bad == 'thermal': ctx['thermal_protocol']['rotor_cutoffs_cm'] = [50]
    if bad == 'graph': ctx['task']['graph_sha256'] = 'other graph'
    if bad == 'component': row['nuclear_by_cutoff_eV']['50'] += .01
    if bad == 'geometry': (Path(row['namespace']).parent/'input.xyz').write_text('changed coordinates')
    if bad == 'receipt': Path(row['review']).write_text('changed receipt')
    if bad == 'roster': native['states'][-1]['key'] = 'invented input'
    if bad == 'audit': audit['producer_sha256'] = 'another producer'
    if bad == 'batches': contract['first_batch_input_IDs'] = [2, 3]
    with pytest.raises(ValueError): run(paths, native, audit, contract)
