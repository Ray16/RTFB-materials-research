"""Collect exact, balanced finite-minimum native occupancy free energies.

Only independently reviewed components with actual SCF-frame checks enter a
cycle. Missing components stay blank; E, G, sampling and material admission have
separate meanings. No new quantum calculation is launched by this module.
"""
from copy import deepcopy
import json
import math
from pathlib import Path

from redox.core.common import ROOT, HARTREE_EV, EV_KJ
from redox.core.hashing import sha256_file
from redox.properties.context import combine_properties
from redox.workflows.ionic_stationary_cycles import validate_cycle_balance, screening_signature
from redox.workflows.molecular_native_namespace import verify, load

CUTOFFS = ('25', '50', '100')


def validate_component(row, expected_charge):
    if (row.get('native_G_qualified') is not True or row.get('conformer_ensemble_qualified')
            or row.get('primary_model_admission') or row.get('material_admission')
            or row['charge'] != expected_charge or row['mult'] != 1
            or set(row['G_by_cutoff_eV']) != set(CUTOFFS)
            or set(row['nuclear_by_cutoff_eV']) != set(CUTOFFS)):
        raise ValueError('Exact finite-minimum component state/scope required')
    for cutoff in CUTOFFS:
        terms = (row['E_geometry_level_Ha'], row['G_by_cutoff_eV'][cutoff], row['nuclear_by_cutoff_eV'][cutoff])
        if (any(type(v) not in (int, float) or not math.isfinite(v) for v in terms)
                or abs(terms[0]*HARTREE_EV + terms[2] - terms[1]) > 1e-8):
            raise ValueError('Finite, independently closed E/nuclear/G component required')


def native_cycle(cycle, records, cases):
    validate_cycle_balance(cycle, cases)
    coefficients = cycle['coefficients']
    if sum(coefficients.values()) != 0:
        raise ValueError('This occupancy cycle must conserve component count')
    missing = [key for key in coefficients if key not in records]
    if missing:
        return dict(status='pending_independently_checked_native_components', missing=missing,
                    delta_delta_E_kJ_mol=None, delta_delta_G_by_cutoff_kJ_mol=None,
                    native_cycle_qualified=False, conformer_ensemble_qualified=False)
    checked = {}
    for key in coefficients:
        row = records[key]
        validate_component(row, cases[key]['charge'])
        checked[key] = row
    electronic = math.fsum(coefficients[k]*checked[k]['E_geometry_level_Ha'] for k in coefficients)*HARTREE_EV
    values, nuclear, assemblies = {}, {}, {}
    for cutoff in CUTOFFS:
        inputs = {}
        for key, row in checked.items():
            context = deepcopy(row['context'])
            context['thermal_definition']['selected_rotor_cutoff_cm'] = int(cutoff)
            inputs[key] = dict(value=row['G_by_cutoff_eV'][cutoff], context=context,
                state=dict(identity=cases[key]['smiles'], charge=row['charge'], multiplicity=row['mult']),
                qc='independently_checked', sources=row['sources'])
        assembly = combine_properties(inputs, coefficients)
        thermal = math.fsum(coefficients[k]*checked[k]['nuclear_by_cutoff_eV'][cutoff] for k in coefficients)
        if abs(assembly['value'] - electronic - thermal) > 1e-8:
            raise ValueError('Balanced E/nuclear/G cycle decomposition does not close')
        values[cutoff] = assembly['value']*EV_KJ
        nuclear[cutoff] = thermal*EV_KJ
        assemblies[cutoff] = assembly
    return dict(status='finite_minimum_native_occupancy_cycle_checked', missing=[],
        delta_delta_E_kJ_mol=electronic*EV_KJ, delta_delta_G_by_cutoff_kJ_mol=values,
        nuclear_contribution_by_cutoff_kJ_mol=nuclear, property_assemblies=assemblies,
        native_cycle_qualified=True, conformer_ensemble_qualified=False,
        primary_model_admission=False, isolated_motif_binding_energy=False,
        electronic_coupling_computed=False, material_admission=False,
        component_count_change=0, one_M_conversion_cancels_for_balanced_cycle=True)


def analyze(preparation, reviews):
    data = json.loads(preparation.read_text())
    auditor = load(ROOT/'scripts/validation/check_charged_polymer_preparation.py', 'native_cycle_exact_preparation')
    auditor.check(preparation)
    cases = {r['case']['key']: r['case'] for r in data['cases']}
    if len(cases) != len(data['cases']):
        raise ValueError('Distinct exact graph/configuration cases required')
    records, protocols = {}, {}
    hashes = {str(preparation): sha256_file(preparation)}
    for path in reviews:
        report = json.loads(path.read_text())
        if (report.get('status') != 'charged_polymer_native_component_independently_reconstructed'
                or report.get('native_G_qualified') is not True
                or type(report.get('actual_SCF_checkpoint_frames_checked')) is not int
                or report['actual_SCF_checkpoint_frames_checked'] < 1
                or report.get('primary_model_admission') or report.get('material_admission')
                or report.get('conformer_ensemble_qualified') or not report.get('source_sha256')):
            raise ValueError('Independent native component with actual checkpoint-frame review required')
        for source, digest in report['source_sha256'].items():
            if sha256_file(source) != digest:
                raise ValueError('Native component review backing changed: '+source)
        namespace = Path(report['namespace']); spec = verify(namespace)
        case, sources = Path(spec['case']), Path(spec['sources'])
        plan = json.loads((case/'plan.json').read_text()); key = report['graph_key']
        if (key in records or key not in cases or plan['prepared_case'] != cases[key]
                or report['candidate_alias'] != data['candidate_alias']
                or report['charge'] != cases[key]['charge'] or report['mult'] != 1
                or report['model'] != plan['solvent_model']
                or report['quantum_settings'] != plan['numerical_settings']
                or plan['calculation_ion_convention'] != data['config']['ionic_model']['calculation_convention']
                or plan['explicit_counterions_in_coordinates'] is not False
                or (not cases[key]['occupied_sites'] and not report.get('matched_empty_polymer_reference'))):
            raise ValueError('Distinct exact graph/state/model/ion conventions required')
        raw = json.loads((case/'thermodynamic_component.json').read_text())
        validate_component(report, cases[key]['charge'])
        if (raw['geometry_sha256'] != report['geometry_sha256']
                or abs(raw['E_geometry_level_Ha']-report['E_geometry_level_Ha'])>1e-10
                or any(abs(raw[field][k]-report[field][k])>1e-8
                       for field in ('G_by_cutoff_eV','nuclear_by_cutoff_eV') for k in CUTOFFS)):
            raise ValueError('Native review differs from the raw thermodynamic component')
        implementation_path = case/'screening_implementation/implementation.json'
        compiled = implementation_path.with_name('screening_kernel.cu')
        signature = screening_signature(json.loads(implementation_path.read_text()), sources, compiled)
        signature['screening_helper_sha256'] = sha256_file(sources/'helper.py')
        qrrho_source = Path(spec['archive'])/'src/redox/qm/qrrho.py'
        digest = sha256_file(qrrho_source)
        if plan['source_sha256'].get(str(qrrho_source)) != digest:
            raise ValueError('Actual native thermal implementation must be source-bound')
        if digest not in protocols:
            protocols[digest] = load(qrrho_source, 'native_cycle_thermal_'+digest[:12]).PROTOCOL
        protocol = deepcopy(protocols[digest])
        if protocol['temperature_K'] != 298.15 or protocol['rotor_cutoffs_cm'] != [25, 50, 100]:
            raise ValueError('Explicit current thermal temperature/cutoffs required')
        context = dict(quantity='native_component_G', unit='eV',
            electronic_protocol=deepcopy(plan['numerical_settings']),
            thermal_definition=dict(protocol=protocol, source_sha256=digest, spin_partition_eV=0.),
            solvent='acetonitrile', solvent_model=deepcopy(plan['solvent_model']),
            temperature_K=protocol['temperature_K'],
            standard_state=dict(translation_convention=protocol['pressure'],
                                one_M_correction_applied=False),
            reference=dict(energy='absolute component within one declared model',
                           ion_convention=plan['calculation_ion_convention'],
                           explicit_spectator_ions=False, numerical_implementation=signature))
        bound = {str(path): sha256_file(path), str(namespace): sha256_file(namespace),
                 str(case/'plan.json'): sha256_file(case/'plan.json')}
        records[key] = {**report, 'context': context, 'sources': bound}
        hashes.update(bound)
        hashes[str(implementation_path)] = sha256_file(implementation_path)
        hashes[str(compiled)] = sha256_file(compiled)
    cycles = [dict(definition=c, **native_cycle(c, records, cases)) for c in data['cycles']]
    count = sum(c['native_cycle_qualified'] for c in cycles)
    return dict(status='exact_ionic_native_occupancy_cycles_collected',
        candidate_alias=data['candidate_alias'], records=records, cycles=cycles,
        checked_components=len(records), computed_cycles=count,
        all_declared_native_cycles_qualified=bool(cycles) and count==len(cycles),
        conformer_ensemble_qualified=False, primary_model_admission=False,
        experimental_loading_inferred=False, material_admission=False,
        source_sha256=hashes, producer_sha256=sha256_file(__file__))


def main():
    import argparse
    import os
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--preparation',type=Path,required=True)
    p.add_argument('--review',type=Path,action='append',required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier native-cycle analyses')
    result=analyze(a.preparation.resolve(),[r.resolve() for r in a.review])
    atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','checked_components','computed_cycles','all_declared_native_cycles_qualified')}))


if __name__=='__main__':main()
