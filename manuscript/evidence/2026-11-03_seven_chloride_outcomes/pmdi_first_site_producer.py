#!/usr/bin/env python
"""Reconstruct available first-graft G and base-canceling polymer penalties.

Missing double-occupancy thermochemistry stays missing. Electronic cycles,
first-site grafting, placement and capped-model differences are distinct outputs.
"""
import argparse
from collections import Counter
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record, record_path
from redox.qm.native_thermal import cache_path, load_native_thermal
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

EV_KJ = 96.48533212331001


def balance(reactants, products):
    from rdkit import Chem
    def inventory(smiles):
        elements, charge = Counter(), 0
        for text in smiles:
            molecule = Chem.AddHs(Chem.MolFromSmiles(text))
            elements.update(a.GetAtomicNum() for a in molecule.GetAtoms())
            charge += Chem.GetFormalCharge(molecule)
        return elements, charge
    left, right = inventory(reactants), inventory(products)
    if left != right:
        raise ValueError('Reaction fails independent atom/charge balance')
    return dict(atom_balanced=True, charge_balanced=True, reactants=reactants, products=products)


def state_values(state, charge, mult, expected_G=None):
    equations = api.load_module('first_graft_independent_nuclear_equations', ROOT / 'scripts/validation/audit/native_partition.py')
    energy, native = load_record(state, charge, mult), load_native_thermal(state, charge, mult)
    surface, _ = cache_path(state, charge, mult)
    independent = equations.read_partition(surface, state / 'opt.xyz')
    if energy is None or native is None or independent is None:
        raise ValueError('Normally admitted and independently reconstructed thermal parents required')
    G = energy['e_smd_eV'] + independent['g_thermal_total_eV']
    if expected_G is not None and (not math.isfinite(expected_G) or abs(G - expected_G) > 1e-8):
        raise ValueError('Stored component G does not reconstruct')
    sources = {str(p): api.digest(p) for p in (record_path(state, charge, mult), surface, state / 'opt.xyz')}
    sources.update(independent['evidence_sha256'])
    return dict(state_dir=str(state), E_eV=energy['e_smd_eV'], G_eV=G,
        G_by_cutoff_eV={str(k): energy['e_smd_eV'] + independent['nuclear_by_cutoff_eV'][str(k)]
            + native['g_spin_partition_eV'] for k in (25, 50, 100)},
        minimum_frequency_cm=min(native['frequencies_cm']), source_sha256=sources)


def deltas(product, empty, capped, local, constants):
    if any(not math.isfinite(v) for v in (product, empty, capped, local, *constants.values())):
        raise ValueError('Finite qualified component G values required')
    first = {name: (product - empty + value) * EV_KJ for name, value in constants.items()}
    cap = {name: (capped - local + value) * EV_KJ for name, value in constants.items()}
    penalty = (product + local - empty - capped) * EV_KJ
    if any(abs(first[name] - cap[name] - penalty) > 1e-7 for name in first):
        raise ValueError('Base-canceling thermodynamic double difference fails closure')
    return dict(first_graft_kJ_mol=first, capped_graft_kJ_mol=cap,
        polymer_minus_capped_kJ_mol=penalty)


def analyze(base, empty_campaign, singles, pair_energy, out):
    from rdkit import Chem
    declaration = json.loads((base / 'baseline_plan.json').read_text())
    api.verify_sources(declaration['source_sha256'])
    snapshot = json.loads((base / 'candidate_snapshot.json').read_text())
    reactions_path = ROOT / 'results/grafting/reactions.json'
    reaction = next(r for r in json.loads(reactions_path.read_text())['reactions'] if r['id'] == snapshot['id'])
    expected_empty = declaration['shared_empty_reference']['graph_smiles']
    empty_rows, invalid = api.admitted_records(empty_campaign, 'thermal')
    matches = [r for r in empty_rows.values() if r['task']['smiles'] == expected_empty]
    if invalid or len(matches) != 2:
        raise ValueError('The declared pair of exact inherited empty references is required')
    empty_record = min(matches, key=lambda r: r['E_smd_eV'])
    empty = state_values(Path(empty_record['state_dir']), 0, 1, empty_record['G_eV'])
    sources = dict(empty['source_sha256'])
    sources.update({str(p): api.digest(p) for p in (Path(__file__).resolve(), base / 'baseline_plan.json',
        base / 'candidate_snapshot.json', reactions_path, empty_campaign / 'plan.json')})
    ref = ROOT / 'results/grafting/campaigns/2026-10-19_mecn_reference'
    base_root = ROOT / 'results/grafting/campaigns/2026-10-19_base_assisted_pilot'
    paths = dict(precursor=ref / 'done' / (reaction['precursor_id'] + '.json'),
        local=ref / 'done/merrifield_local_chloride.json', chloride=ref / 'done/chloride.json',
        base=base_root / 'done/dbu_neutral.json', protonated=base_root / 'done/dbu_protonated.json',
        pair=base_root / 'done/dbu_chloride_pair.json')
    reagents = {}; smiles = {}
    for name, path in paths.items():
        record = verify_state_checkpoint(path); task = record['task']; row = task.get('manifest_row', task)
        text = row['smiles']; molecule = Chem.MolFromSmiles(text)
        q = Chem.GetFormalCharge(molecule); m = int(row.get('mult_hint', task.get('mult', 1)))
        if name == 'chloride':
            if record.get('readiness') != 'qualified_monatomic_reference':
                raise ValueError('Qualified monatomic chloride required')
            value = dict(G_eV=record['G_eV'], G_by_cutoff_eV={str(k): record['G_eV'] for k in (25, 50, 100)},
                state_dir=record['state_dir'], source_sha256=record['output_sha256'])
        else:
            value = state_values(Path(record['state_dir']), q, m, record['G_eV'])
        reagents[name] = value; smiles[name] = text
        sources[str(path)] = api.digest(path); sources.update(record['output_sha256']); sources.update(value['source_sha256'])
    reused_path = ref / 'reused_products.json'; reused = json.loads(reused_path.read_text())[snapshot['id']]
    api.verify_sources(reused['source_sha256'])
    capped = state_values(Path(reused['state_dir']), 0, 1, reused['G_eV'])
    sources[str(reused_path)] = api.digest(reused_path); sources.update(reused['source_sha256']); sources.update(capped['source_sha256'])
    correction = -.00831446261815324 * 298.15 * math.log(.08205736608096 * 298.15)
    constants = {}
    for cutoff in (25, 50, 100):
        g = {name: value['G_by_cutoff_eV'][str(cutoff)] for name, value in reagents.items()}
        constants[str(cutoff)] = dict(free_ions=g['protonated'] + g['chloride'] - g['base'] - g['precursor'],
            contact_pair=g['pair'] - g['base'] - g['precursor'] + correction / EV_KJ)
    old_terms_path = base / 'reaction_terms.json'; old_terms = json.loads(old_terms_path.read_text())
    api.verify_sources(old_terms['source_sha256']); sources[str(old_terms_path)] = api.digest(old_terms_path)
    expected_pair_shift = correction * (1 / EV_KJ - 1 / 96.485)
    if (abs(constants['50']['free_ions'] - old_terms['constant_eV']['free_ions']) > 1e-8
            or abs(constants['50']['contact_pair'] - old_terms['constant_eV']['contact_pair'] - expected_pair_shift) > 1e-8):
        raise ValueError('Shared reagent term differs beyond the explicit conversion-constant update')
    rows = []; lookup = {}
    geometry_audit = api.load_module('first_graft_independent_polymer_geometry', ROOT / 'scripts/validation/audit_polymer_geometries.py')
    for campaign in singles:
        admitted, invalid = api.admitted_records(campaign, 'thermal')
        if invalid or len(admitted) != 1:
            raise ValueError('One normally admitted native component per single-site campaign required')
        record = next(iter(admitted.values())); task = record['task']; occupied = task['occupied_sites']
        if len(occupied) != 1 or occupied[0] in lookup:
            raise ValueError('Distinct single-site occupancy components required')
        value = state_values(Path(record['state_dir']), 0, 1, record['G_eV'])
        audit = geometry_audit.audit(task, Path(record['state_dir']) / 'opt.xyz')
        if audit['issues'] or task['backbone_stereo_fingerprint'] != empty_record['task']['backbone_stereo_fingerprint']:
            raise ValueError('Exact polymer geometry or backbone configuration differs')
        free_balance = balance([expected_empty, smiles['precursor'], smiles['base']],
            [task['smiles'], smiles['protonated'], smiles['chloride']])
        pair_balance = balance([expected_empty, smiles['precursor'], smiles['base']], [task['smiles'], smiles['pair']])
        relative_balance = balance([task['smiles'], reaction['site_smiles']], [expected_empty, reaction['product_smiles']])
        by_cutoff = {str(k): deltas(value['G_by_cutoff_eV'][str(k)], empty['G_by_cutoff_eV'][str(k)],
            capped['G_by_cutoff_eV'][str(k)], reagents['local']['G_by_cutoff_eV'][str(k)], constants[str(k)]) for k in (25, 50, 100)}
        row = dict(site=occupied[0], task=task, state_dir=value['state_dir'], E_eV=value['E_eV'], G_eV=value['G_eV'],
            state_G_by_cutoff_eV=value['G_by_cutoff_eV'],
            minimum_frequency_cm=value['minimum_frequency_cm'], by_cutoff=by_cutoff,
            balances=dict(free_ions=free_balance, contact_pair=pair_balance, polymer_minus_capped=relative_balance),
            geometry_audit=audit)
        rows.append(row); lookup[occupied[0]] = row
        sources.update(value['source_sha256']); sources.update(record['output_sha256'])
        sources[str(campaign / 'plan.json')] = api.digest(campaign / 'plan.json')
    if set(lookup) != {2, 3}:
        raise ValueError('This declared comparison requires sites 2 and 3')
    paired, invalid = api.admitted_records(pair_energy, 'energy')
    if invalid or len(paired) != 1:
        raise ValueError('The paired electronic record must pass ordinary admission')
    pair = next(iter(paired.values()))
    if pair['task']['occupied_sites'] != [2, 3]:
        raise ValueError('Wrong double-occupancy electronic model')
    balance([pair['task']['smiles'], expected_empty], [lookup[2]['task']['smiles'], lookup[3]['task']['smiles']])
    sources.update(pair['output_sha256']); sources[str(pair_energy / 'plan.json')] = api.digest(pair_energy / 'plan.json')
    placement = dict(site2_minus_site3_deltaE_kJ_mol=(lookup[2]['E_eV'] - lookup[3]['E_eV']) * EV_KJ,
        site2_minus_site3_deltaG_kJ_mol=(lookup[2]['G_eV'] - lookup[3]['G_eV']) * EV_KJ)
    report = dict(status='available_first_site_grafting_independently_reconstructed', checked=time.time(),
        candidate_id=snapshot['id'], solvent='acetonitrile', temperature_K=298.15, base='DBU',
        conversion_kJ_mol_per_eV=EV_KJ, contact_pair_standard_state_correction_kJ_mol=correction,
        shared_constant_eV_by_cutoff=constants, singles=sorted(rows, key=lambda r: r['site']),
        reagent_states={name: dict(state_dir=value['state_dir'], G_eV=value['G_eV'],
            G_by_cutoff_eV=value['G_by_cutoff_eV'], smiles=smiles[name]) for name, value in reagents.items()},
        empty_reference=dict(state_dir=empty['state_dir'], E_eV=empty['E_eV'], G_eV=empty['G_eV'],
            G_by_cutoff_eV=empty['G_by_cutoff_eV'], smiles=expected_empty,
            selection='Lowest electronic-energy member of the two inherited exact-graph empty conformers'),
        local_capped_reference=dict(state_dir=capped['state_dir'], G_eV=capped['G_eV'],
            G_by_cutoff_eV=capped['G_by_cutoff_eV'], smiles=reaction['product_smiles']),
        placement=placement, electronic_neighbor_cycle_kJ_mol=(pair['E_smd_eV'] + empty['E_eV'] - lookup[2]['E_eV'] - lookup[3]['E_eV']) * EV_KJ,
        paired_G_eV=None, neighbor_deltaDeltaG_kJ_mol=None, conformer_convergence_established=False,
        synthesis_admission=False, source_sha256=sources,
        scope='Finite one-minimum first-site reactions and common-reagent-canceling local/polymer differences; no complete neighbor free-energy cycle, loading, isolated motif binding or coupling')
    api.verify_sources(sources)
    out.mkdir(parents=True, exist_ok=False); atomic_write_json(out / 'summary.json', report)
    print(json.dumps(dict(placement=placement, first_sites=[dict(site=r['site'], **r['by_cutoff']['50']) for r in report['singles']])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, required=True)
    parser.add_argument('--empty-campaign', type=Path, required=True)
    parser.add_argument('--single-campaign', type=Path, action='append', required=True)
    parser.add_argument('--pair-energy', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    analyze(args.base.resolve(), args.empty_campaign.resolve(), [p.resolve() for p in args.single_campaign],
        args.pair_energy.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
