"""Analyze complete corrected native cells after a source-bound snapshot audit.

Retains all twelve declared cells. Missing inputs leave that cell blank; a
complete cell does not qualify a balanced cycle or an experimental ensemble.
"""
import argparse
import importlib.metadata
import json
import math
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.properties.conformer_geometry import alignment_molecule, symmetry_rmsd_matrix
from redox.validation.native_conformer_batches import finite_batches
from redox.validation.sampling_snapshot import roster

COMMON_CONTEXT = ('experimental_method_id', 'experimental_model', 'electronic_protocol',
                  'native_settings', 'thermal_protocol', 'standard_state', 'temperature_K',
                  'charge', 'multiplicity')
GRAPH_CONTEXT = ('smiles', 'mapped_smiles', 'graph_sha256', 'occupied_sites', 'reactive_sites',
                 'site_states', 'backbone_stereo_fingerprint', 'backbone_carbon_count',
                 'backbone_CC_bonds', 'chain_ends', 'atoms', 'charge', 'mult', 'model_key')


def common_context(row, contract):
    context = row['context']
    required = {k: contract[k] for k in ('experimental_method_id', 'electronic_protocol', 'native_settings')}
    required.update(experimental_model=contract['solvent_model'], charge=0, multiplicity=1,
                    temperature_K=298.15)
    if (any(context[k] != value for k, value in required.items())
            or context['thermal_protocol']['temperature_K'] != context['temperature_K']
            or context['thermal_protocol']['rotor_cutoffs_cm'] != contract['rotor_cutoffs_cm']
            or context['thermal_protocol']['pressure'] != context['standard_state']):
        raise ValueError('Exact declared native, electronic, solvent, state and thermal context required')
    return {k: context[k] for k in COMMON_CONTEXT}


def checked_components(row, audited):
    if (row['status'] != 'native_independently_checked' or row['native_G_qualified'] is not True
            or row['G_by_cutoff_eV'] != audited['G_by_cutoff_eV']
            or row['review'] != audited['review']):
        raise ValueError('Native component differs from the audited state')
    if set(row['G_by_cutoff_eV']) != {'25', '50', '100'} or set(row['nuclear_by_cutoff_eV']) != {'25', '50', '100'}:
        raise ValueError('All three declared rotor settings required')
    electronic = row['electronic_E_smd_eV']
    for cutoff, free in row['G_by_cutoff_eV'].items():
        nuclear = row['nuclear_by_cutoff_eV'][cutoff]
        if (any(type(v) not in (float, int) or not math.isfinite(v) for v in (electronic, nuclear, free))
                or abs(free - electronic - nuclear) > 1e-8):
            raise ValueError('Native G must retain its own finite checked electronic and nuclear terms')


def analyze(native_path, audit_path, contract_path, *, source_reader=None):
    def verify_sources(values):
        for path,digest in values.items():
            if source_reader is not None:
                source_reader.verify(path,digest)
            elif sha256_file(path)!=digest:
                raise ValueError('An audited input or receipt changed')
    if source_reader is not None:
        source_reader.allow_local((native_path,audit_path,contract_path))
    sources = {str(p): sha256_file(p) for p in (native_path, audit_path, contract_path)}
    audit = json.loads(audit_path.read_text()); native = json.loads(native_path.read_text())
    contract = json.loads(contract_path.read_text())
    rows = roster(native, 'complete_roster_corrected_native_review_snapshot')
    if (audit['status'] != 'frozen_sampling_component_receipts_rechecked'
            or audit['producer_sha256'] != sha256_file(Path(__file__).with_name('sampling_snapshot.py'))
            or sources[str(native_path)] not in audit['source_sha256'].values()
            or not audit['context_source_and_component_sums_rechecked']
            or audit['occupied_ensemble_cycles_qualified'] or audit['primary_candidate_admission']
            or contract['global_convergence_proven'] or contract['historical_E_or_G_imported']):
        raise ValueError('A matching current source-bound component audit and unqualified sampling contract required')
    sources.update(audit['source_sha256'])
    verify_sources(sources)
    indices = contract['expanded_batch_input_IDs']; early = contract['first_batch_input_IDs']
    if (indices != list(range(1, contract['starts_per_cell'] + 1)) or early != list(range(1, len(early) + 1))
            or not 0 < len(early) < len(indices)):
        raise ValueError('Contiguous strictly nested declared batches required')
    cells = {}
    for configuration in contract['configurations']:
        for sites in contract['occupancy_graphs']:
            key = configuration + '__occ_' + ('_'.join(map(str, sites)) if sites else 'none')
            cells[key] = [key + f'__conf{i:03d}' for i in indices]
    if set(sum(cells.values(), [])) != set(rows):
        raise ValueError('Snapshot must retain the exact full declared graph and input roster')
    audited = {r['key']: r for r in audit['records']}
    checked = {k: r for k, r in rows.items() if r['native_G_qualified'] is True}
    if len(audited) != len(audit['records']) or set(audited) != set(checked):
        raise ValueError('Audit and snapshot native qualification differ')
    if len(checked) != native['checked_native_inputs'] or len(checked) != audit['checked_native_inputs']:
        raise ValueError('Checked native census differs')
    global_context = None
    for key, row in checked.items():
        checked_components(row, audited[key])
        context = common_context(row, contract)
        if global_context is not None and context != global_context:
            raise ValueError('Different physical contexts cannot be combined')
        global_context = context
    results = []; bindings = {}
    for key, expected in cells.items():
        missing = [k for k in expected if k not in checked]
        cell = dict(model_key=key, checked_inputs=len(expected)-len(missing), missing_keys=missing)
        if missing:
            cell.update(status='pending_complete_native_batch', partitions=None)
        else:
            selected = [checked[k] for k in expected]; molecules = []; graph = None
            for row in selected:
                task = row['context']['task']; identity = {k: task[k] for k in GRAPH_CONTEXT}
                if (task['key'] != row['key'] or task['model_key'] != key
                        or (graph is not None and identity != graph)):
                    raise ValueError('Exact graph, backbone and site identities must agree within each cell')
                graph = identity
                namespace_path = Path(row['namespace'])
                namespace_input = source_reader.path(namespace_path) if source_reader is not None else namespace_path
                sources[str(namespace_path)] = sha256_file(namespace_input)
                spec = json.loads(namespace_input.read_text()); geometry = Path(spec['case'])/'input.xyz'
                sources[str(geometry)] = row['context']['geometry_sha256']
                geometry_input = source_reader.verify(geometry,sources[str(geometry)]) if source_reader is not None else geometry
                if sha256_file(geometry_input) != sources[str(geometry)]:
                    raise ValueError('Native geometry differs from the checked electronic/nuclear input')
                molecules.append(alignment_molecule(task['smiles'], geometry_input))
                bindings[row['key']] = dict(namespace=str(namespace_path), geometry=str(geometry),
                    geometry_sha256=sources[str(geometry)], review=row['review'], review_sha256=row['review_sha256'])
            matrix, matches = symmetry_rmsd_matrix(molecules)
            records = [dict(key=r['key'], E_eV=r['electronic_E_smd_eV'],
                            G_by_cutoff_eV=r['G_by_cutoff_eV']) for r in selected]
            cell.update(status='complete_finite_native_batches_compared', graph_context=graph,
                rmsd_A=matrix.tolist(), symmetry_matches=matches,
                partitions=finite_batches(records, matrix, len(early), temperature_K=global_context['temperature_K'],
                    rmsd_cutoffs_A=contract['duplicate_heavy_atom_RMSD_cutoffs_A'],
                    energy_tolerance_kJ_mol=contract['duplicate_electronic_energy_tolerance_kJ_mol'],
                    rotor_cutoffs_cm=contract['rotor_cutoffs_cm'],
                    tolerance_kJ_mol=contract['finite_batch_change_tolerance_kJ_mol']))
        results.append(cell)
    for module in ('sampling_batches.py', 'native_conformer_batches.py', 'sampling_snapshot.py'):
        p = Path(__file__).with_name(module); sources[str(p)] = sha256_file(p)
        if source_reader is not None:source_reader.allow_local((p,))
    from redox.properties import conformer_geometry, conformer_ensemble
    from redox.build import geometry as geometry_code
    from redox.validation import protonation
    for module in (conformer_geometry, conformer_ensemble, geometry_code, protonation):
        p = Path(module.__file__); sources[str(p)] = sha256_file(p)
        if source_reader is not None:source_reader.allow_local((p,))
    verify_sources(sources)
    return dict(status='audited_corrected_native_finite_batches_analyzed', declared_inputs=len(rows),
        checked_native_inputs=len(checked), complete_cells=sum(r['partitions'] is not None for r in results),
        first_count=len(early), total_count=len(indices), cells=results, context=global_context, inputs=bindings,
        duplicate_policy='complete-link full-set numerical well assignment restricted to nested batches; minimum-G representative',
        duplicate_energy='each native component own checked wb97m-v/def2-tzvpd SMD electronic parent',
        alignment_policy='graph and hydrogen checked; heavy-atom symmetry maps preserve chirality; proper rotations only',
        versions={p: importlib.metadata.version(p) for p in ('numpy', 'rdkit', 'ase')},
        original_physical_QC_repeated=False, quantum_calculations_run=False,
        occupied_ensemble_cycles_qualified=False, global_ensemble_convergence_proven=False,
        experimental_occupancy_distribution_assumed=False, primary_candidate_admission=False,
        material_admission=False, source_sha256=sources)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'audit', 'contract', 'out'): p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    if a.out.exists(): raise FileExistsError('Preserve prior finite-batch analyses')
    result = analyze(a.native.resolve(), a.audit.resolve(), a.contract.resolve())
    atomic_write_json(a.out, result)
    print(json.dumps({k: result[k] for k in ('status', 'checked_native_inputs', 'complete_cells')}))
