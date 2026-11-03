#!/usr/bin/env python3
"""Assemble the exact refined phenothiazine/Fc cycle from qualified components."""
import argparse
import ast
from copy import deepcopy
from decimal import Decimal
import json
import math
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.runtime_paths import gpu4pyscf_source_root
from redox.properties.reference_redox import COEFFICIENTS, CUTOFFS, neutral_oxidation_vs_fc
from redox.workflows.molecular_native_namespace import verify


def assemble(output):
    from ase.io import read as read_atoms
    if output.exists():
        raise FileExistsError('Preserve earlier redox assemblies')
    used = {}
    def read(path):
        path = Path(path); used[str(path)] = sha256_file(path)
        return json.loads(path.read_text())
    def checked_review(path):
        review = read(path)
        for name, digest in review['source_sha256'].items():
            if sha256_file(name) != digest:
                raise ValueError('Changed component backing: ' + name)
            used[name] = digest
        return review

    folder = ROOT / 'validation/phenothiazine_cds_2026-11-03'
    contract = read(ROOT / 'validation/phenothiazine_primary_2026-11-02/calculation_contract.json')
    reference_path = Path(contract['reference']['file'])
    if sha256_file(reference_path) != contract['reference']['sha256']:
        raise ValueError('Declared Fc reference cohort changed')
    primary = read(reference_path)
    contexts = [primary['components'][key]['context'] for key in ('fc_neu', 'fc_ox')]
    if contexts[0] != contexts[1]:
        raise ValueError('Fc reference contexts differ')
    context = deepcopy(contexts[0])
    ptz = {key: checked_review(folder / leaf) for key, leaf in (
        ('minus', 'neutral_minus_native_review.json'),
        ('plus', 'neutral_plus_native_review.json'),
        ('oxidized', 'oxidized_finer_native_review.json'))}
    thermal = ptz['minus']['thermal_protocol']
    # Retain the complete kernel/mass/projection/standard-state specification.
    context['thermal_definition']['full_projected_protocol'] = thermal
    cds_source = gpu4pyscf_source_root() / 'solvent/hessian/smd.py'
    tree = ast.parse(cds_source.read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'get_cds')
    steps = [ast.literal_eval(n.value) for n in function.body if isinstance(n, ast.Assign)
             and any(isinstance(t, ast.Name) and t.id == 'eps' for t in n.targets)]
    if steps != [1e-4]:
        raise ValueError('Validated default CDS implementation changed')
    used[str(cds_source)] = sha256_file(cds_source)
    default_policy = dict(schema=1, component='primary legacy SMD CDS',
        Hessian_gradient_step_Bohr=steps[0], energy_function_changed=False,
        gradient_function_changed=False, global_defaults_changed=False)
    records = {}

    def record(key, review, review_path, policy):
        spec = verify(Path(review['namespace'])); case = Path(spec['case'])
        plan = read(case / 'plan.json')
        energy = read(Path(spec['sources']) / 'corrected_electronic.json')
        native_paths = list((case / 'state/native_thermal').glob('*.json'))
        projection_paths = list((case / 'state/projected_thermal').glob('*.json'))
        if len(native_paths) != 1 or len(projection_paths) != 1:
            raise ValueError('One own native and projected record required')
        native = read(native_paths[0]); projection = read(projection_paths[0])
        if (energy['energy_protocol'] != context['electronic_protocol']
                or native['identity']['protocol'] != context['thermal_definition']['native_protocol']
                or projection['identity']['protocol'] != thermal
                or projection['kernel_execution']['library_sha256'] != thermal['library_sha256']
                or review['numerical_implementation'] != context['reference']['numerical_implementation']
                or (review['charge'], review['mult']) != (plan['charge'], plan['mult'])):
            raise ValueError('Actual electronic/native/partition/state contexts differ')
        for cutoff in CUTOFFS:
            value = energy['e_smd_eV'] + projection['nuclear_by_cutoff_eV'][cutoff] + projection['g_spin_partition_eV']
            if not math.isfinite(value) or abs(value-review['G_by_cutoff_eV'][cutoff]) > 1e-8:
                raise ValueError('Own component sum does not reproduce reviewed G')
        geometry = case / 'state/opt.xyz'; used[str(geometry)] = sha256_file(geometry)
        return dict(state=dict(identity=key, graph_identity=plan['smiles'],
                    charge=review['charge'], multiplicity=review['mult']),
            symbols=read_atoms(geometry).get_chemical_symbols(),
            context=deepcopy(context), G_by_cutoff_eV=review['G_by_cutoff_eV'],
            native_numerical_policy=deepcopy(policy), native_G_qualified=True,
            thermodynamic_entities=1, qc='independently_checked',
            sources={str(review_path): sha256_file(review_path)},
            actual_native_SCF_frames_checked=review['actual_SCF_checkpoint_frames_checked'])

    for key in ('fc_neu', 'fc_ox'):
        component = primary['components'][key]
        if len(component['sources']) != 1:
            raise ValueError('One qualified Fc review required')
        path, digest = next(iter(component['sources'].items()))
        if sha256_file(path) != digest:
            raise ValueError('Fc review differs from the declared cohort')
        review = checked_review(path)
        if (review['status'] != 'primary_molecular_native_control_independently_reconstructed'
                or not review['native_component_requalified']
                or review['actual_SCF_checkpoint_frames_checked'] != 9
                or review['actual_SCF_wavefunction_spins_checked'] != 9
                or review['G_by_cutoff_eV'] != component['G_by_cutoff_eV']):
            raise ValueError('Qualified own Fc thermochemistry and actual spins required')
        records[key] = record(key, review, Path(path), default_policy)
    for key, review in ptz.items():
        if (not review['native_G_qualified'] or review['actual_SCF_checkpoint_frames_checked'] != 21
                or review['smiles'] != contract['product_smiles']
                or review['native_numerical_policy']['Hessian_gradient_step_Bohr'] != 1e-6):
            raise ValueError('Exact qualified finer-step phenothiazine states required')
        name = 'oxidized_finer' if key == 'oxidized' else 'neutral_' + key
        records[key] = record(key, review, folder / (name + '_native_review.json'), review['native_numerical_policy'])
    cycles = {}; discrepancies = []
    for branch in ('minus', 'plus'):
        components = dict(graft_neu=records[branch], graft_ox=records['oxidized'],
                          fc_neu=records['fc_neu'], fc_ox=records['fc_ox'])
        cycles[branch] = neutral_oxidation_vs_fc(components)
        for cutoff in CUTOFFS:
            independent = sum(Decimal(c) * Decimal(str(components[k]['G_by_cutoff_eV'][cutoff]))
                              for k, c in COEFFICIENTS.items())
            discrepancies.append(abs(float(independent)-cycles[branch]['potential_by_cutoff_V'][cutoff]))
    if max(discrepancies) > 1e-9:
        raise ValueError('Independent decimal cycle does not close')
    used[str(Path(__file__))] = sha256_file(__file__)
    used[str(Path(neutral_oxidation_vs_fc.__code__.co_filename))] = sha256_file(neutral_oxidation_vs_fc.__code__.co_filename)
    result = dict(status='phenothiazine_finite_oxidation_potential_context_and_arithmetic_checked',
        cycles=cycles, primary_rotor_cutoff_cm=50, neutral_branches_not_population_weighted=True,
        maximum_independent_decimal_error_V=max(discrepancies),
        numerical_policy_note='Qualified 1e-6-Bohr CDS differentiation for phenothiazine; validated unmodified 1e-4-Bohr default for Fc. Same physical model; component policies retained. No Fc step-convergence study inferred.',
        redox_potential_qualified=True, conformer_ensemble_qualified=False,
        experimental_accuracy_validated=False, material_admission=False, source_sha256=used)
    output.write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True); args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    result = assemble(args.out.resolve())
    print(json.dumps({key: value['potential_by_cutoff_V'] for key, value in result['cycles'].items()}))
