"""Balanced finite-pose electronic cycles; never substitute E for missing G."""
from collections import Counter
import json
import math
from pathlib import Path

from rdkit import Chem

from redox.core.common import EV_KJ
from redox.core.hashing import sha256_file
from redox.workflows.molecular_native_namespace import verify


def electronic_cycle(cycle,records,cases):
    elements=Counter();charge=0
    for key,coefficient in cycle['coefficients'].items():
        if type(coefficient) is not int:raise ValueError('Exact integer cycle coefficients required')
        case=cases[key];mol=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        if Chem.GetFormalCharge(mol)!=case['charge']:raise ValueError('Declared formal state differs from graph')
        for atom in mol.GetAtoms():elements[atom.GetSymbol()]+=coefficient
        charge+=coefficient*case['charge']
    if charge or any(elements.values()):raise ValueError('Electronic cycle is not atom/charge balanced')
    missing=[key for key in cycle['coefficients'] if key not in records]
    if missing:return dict(status='pending_checked_stationary_components',missing=missing,delta_delta_E_Ha=None,delta_delta_E_kJ_mol=None,delta_delta_G_kJ_mol=None)
    delta=math.fsum(coefficient*records[key]['E_Ha'] for key,coefficient in cycle['coefficients'].items())
    if not math.isfinite(delta):raise ValueError('Finite electronic cycle required')
    return dict(status='finite_stationary_electronic_cycle_computed',missing=[],delta_delta_E_Ha=delta,
        delta_delta_E_kJ_mol=delta*27.211386245988*EV_KJ,delta_delta_G_kJ_mol=None,
        conformer_ensemble_qualified=False,isolated_motif_binding_energy=False,electronic_coupling_computed=False)


def screening_signature(receipt,sources,compiled_source):
    if (receipt['status']!='run_local_screening_context_closed_qualification_separate'
            or receipt['upstream_fix']!='https://github.com/pyscf/gpu4pyscf/pull/925'
            or receipt['source_sha256']!=sha256_file(sources/'kernel.cu')
            or receipt['kernel_sha256']!=sha256_file(compiled_source)
            or receipt['installed_library_files_modified']):raise ValueError('Exact corrected screening implementation required')
    return {key:receipt[key] for key in ('upstream_fix','kernel_sha256','source_sha256','installed_gen_grid_sha256')}


def analyze(preparation,reviews):
    data=json.loads(preparation.read_text());cases={r['case']['key']:r['case'] for r in data['cases']};records={};model=None;settings=None;implementation=None;hashes={str(preparation):sha256_file(preparation)}
    for path in reviews:
        report=json.loads(path.read_text())
        if report['status']!='charged_polymer_solution_stationary_geometry_independently_checked':raise ValueError('Independent stationary structure review required')
        for source,digest in report['source_sha256'].items():
            if sha256_file(source)!=digest:raise ValueError('Stationary cycle input changed: '+source)
        namespace=Path(report['namespace']);spec=verify(namespace);plan=json.loads((Path(spec['case'])/'plan.json').read_text());key=report['graph_key']
        if (key in records or key not in cases or plan['prepared_case']!=cases[key]
                or report['candidate_alias']!=data['candidate_alias'] or report['charge']!=cases[key]['charge'] or report['mult']!=1):
            raise ValueError('Distinct exact matched graph/state/candidate inputs required')
        if (plan['calculation_ion_convention']!=data['config']['ionic_model']['calculation_convention']
                or plan['explicit_counterions_in_coordinates'] is not False):raise ValueError('Matched ionic conventions required')
        if not cases[key]['occupied_sites'] and not plan.get('matched_empty_polymer_reference'):
            raise ValueError('An explicit matched empty-reference calculation is required')
        if model is None:model=plan['solvent_model'];settings=plan['quantum_settings']
        if plan['solvent_model']!=model or plan['quantum_settings']!=settings:raise ValueError('Different numerical models cannot form one electronic cycle')
        implementation_path=Path(spec['case'])/'screening_implementation/implementation.json'
        compiled_source=implementation_path.with_name('screening_kernel.cu')
        signature=screening_signature(json.loads(implementation_path.read_text()),Path(spec['sources']),compiled_source)
        if implementation is None:implementation=signature
        if signature!=implementation:raise ValueError('Different screening implementations cannot form one cycle')
        if not math.isfinite(report['E_Ha']):raise ValueError('Finite checked geometry-level electronic energy required')
        records[key]=dict(E_Ha=report['E_Ha'],charge=report['charge'],geometry_sha256=report['geometry_sha256'],review=str(path),review_sha256=sha256_file(path),native_G_qualified=False)
        hashes[str(path)]=sha256_file(path);hashes[str(namespace)]=sha256_file(namespace)
        hashes[str(implementation_path)]=sha256_file(implementation_path)
        hashes[str(compiled_source)]=sha256_file(compiled_source)
    cycles=[dict(definition=cycle,**electronic_cycle(cycle,records,cases)) for cycle in data['cycles']]
    return dict(status='exact_ionic_stationary_electronic_cycles_collected',candidate_alias=data['candidate_alias'],model=model,
        quantum_settings=settings,numerical_implementation=implementation,records=records,cycles=cycles,computed_cycles=sum(row['status']=='finite_stationary_electronic_cycle_computed' for row in cycles),
        native_G_qualified=False,conformer_ensemble_qualified=False,experimental_loading_inferred=False,material_admission=False,
        source_sha256=hashes,producer_sha256=sha256_file(__file__))
