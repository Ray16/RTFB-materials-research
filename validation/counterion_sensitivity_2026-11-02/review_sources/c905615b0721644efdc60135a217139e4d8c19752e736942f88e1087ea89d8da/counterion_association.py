"""Read verified microstates and matched primary parents without model mixing."""
from copy import deepcopy
import json
from pathlib import Path

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.properties.counterion import association_cycle, COEFFICIENTS
from redox.workflows.molecular_native_namespace import verify

SOURCE_BYTES = Path(__file__).read_bytes()


def analyze(primary_inputs, geometry_manifest, component_reviews, output):
    if output.exists():
        raise FileExistsError('Preserve earlier association snapshots')
    data = json.loads(primary_inputs.read_text())
    manifest = json.loads(geometry_manifest.read_text())
    cohort = data['package']
    originals = {s['id']:s for s in cohort['states']}
    context = deepcopy(data['components']['graft_ox2']['context'])
    if context != data['components']['chloride']['context']:
        raise ValueError('Isolated graft and chloride source contexts differ')
    from redox.validation.execution_review import freeze_source
    producers = [freeze_source(output.parent, Path(__file__), SOURCE_BYTES),
        freeze_source(output.parent, ROOT/'src/redox/properties/counterion.py'),
        freeze_source(output.parent, ROOT/'src/redox/properties/context.py')]
    sources = {str(p):sha256_file(p) for p in (primary_inputs, geometry_manifest, *producers)}
    parent_records = {}
    for source_id, role in (('graft_ox2','graft'), ('chloride','chloride')):
        component = data['components'][source_id]
        if component['qc'] != 'independently_checked' or not component['native_component_requalified']:
            raise ValueError('Fully reviewed primary parent components required')
        scalar_backing = []
        for path, digest in component['sources'].items():
            if sha256_file(path) != digest:
                raise ValueError('Primary parent component source changed')
            sources[path] = digest
            backing = json.loads(Path(path).read_text())
            if backing.get('status') in ('primary_molecular_native_control_independently_reconstructed',
                                          'atomic_translation_and_saved_wavefunction_independently_checked'):
                scalar_backing.append(backing)
            for nested, value in backing.get('source_sha256', {}).items():
                if sha256_file(nested) != value:
                    raise ValueError('Primary parent physical review backing changed')
                sources[nested] = value
        if (len(scalar_backing) != 1 or scalar_backing[0]['G_by_cutoff_eV'] != component['G_by_cutoff_eV']
                or (scalar_backing[0]['charge'], scalar_backing[0]['mult']) != (component['charge'], component['mult'])):
            raise ValueError('Primary scalar values or state differ from their physical review')
        parent_records[role] = dict(G_by_cutoff_eV=component['G_by_cutoff_eV'], context=deepcopy(context),
            state=dict(identity=source_id, smiles=originals[source_id]['smiles'], charge=component['charge'],
                       multiplicity=component['mult']), qc='independently_checked', sources=component['sources'],
            thermodynamic_entities=1)
    # The scalar parent convention describes separated ions. For this explicit
    # cycle, state identities distinguish the free ions from the single cluster;
    # the common reference labels the same absolute energies and omitted I-.
    reference = dict(energy=context['reference']['energy'],
        ions='explicit incoming chloride; omitted separated iodide spectator on both sides',
        numerical_implementation=context['reference']['numerical_implementation'])
    for record in parent_records.values():
        record['context']['reference'] = deepcopy(reference)
    rows = {r['key']:dict(key=r['key'], status='pending_checked_microstate_thermochemistry',
            source_convention_deltaG_kJ_mol=None, one_M_deltaG_kJ_mol=None) for r in manifest['states']}
    seen = set()
    for path in component_reviews:
        from redox.workflows.counterion_thermochemistry import checked_component
        report, spec = checked_component(path); key = report['key']
        if key not in rows or key in seen:
            raise ValueError('Distinct declared microstate starts required')
        seen.add(key)
        case = Path(spec['case']); plan = json.loads((case/'plan.json').read_text()); raw = json.loads((case/'result.json').read_text())
        canonical = lambda value: Chem.MolToSmiles(Chem.MolFromSmiles(value))
        if (report['electronic_protocol'] != context['electronic_protocol']
                or report['native_protocol'] != context['thermal_definition']['native_protocol']
                or report['temperature_K'] != context['temperature_K']
                or report['numerical_implementation'] != context['reference']['numerical_implementation']
                or report['standard_state'] != report['thermal_protocol']['pressure']
                or canonical(report['smiles']) != canonical(originals['graft_ox2']['smiles']+'.[Cl-]')
                or (report['charge'], report['mult']) != (1, 1)
                or plan['geometry_sha256'] != report['geometry_sha256']
                or any(abs(value - raw['G_by_cutoff_eV'][c]) > 1e-8 for c, value in report['G_by_cutoff_eV'].items())):
            raise ValueError('Pair and separated parents have incompatible graph, state, model or raw G')
        # Require the full projected thermal protocol used by the actual parent,
        # including the pinned mass/kernel/translation definitions.
        parent_review_path = next(iter(data['components']['graft_ox2']['sources']))
        parent_review = json.loads(Path(parent_review_path).read_text())
        parent_spec = verify(Path(parent_review['namespace']))
        parent_case = Path(parent_spec['case'])
        parent_result = json.loads((parent_case/'result.json').read_text())
        native_path = Path(parent_result['native_record_path'])
        partitions = [p for p in (parent_case/'state/projected_thermal').glob('*.json')
            if json.loads(p.read_text()).get('identity', {}).get('raw_native_sha256') == sha256_file(native_path)]
        if len(partitions) != 1 or json.loads(partitions[0].read_text())['identity']['protocol'] != report['thermal_protocol']:
            raise ValueError('Actual parent and pair projected nuclear protocols differ')
        pair_context = deepcopy(context); pair_context['reference'] = deepcopy(reference)
        paired = dict(G_by_cutoff_eV=report['G_by_cutoff_eV'], context=pair_context,
            state=dict(identity='oxidized_graft_chloride:'+key, smiles=report['smiles'], charge=1, multiplicity=1),
            qc='independently_checked', sources={str(path):sha256_file(path)}, thermodynamic_entities=1)
        cycle = association_cycle(dict(paired=paired, **parent_records))
        rows[key] = dict(key=key, **cycle)
        sources[str(path)] = sha256_file(path); sources.update(report['source_sha256'])
    result = dict(status='primary_counterion_association_review_snapshot', declared_microstate_starts=len(rows),
        computed_microstate_cycles=len(seen), states=list(rows.values()), coefficients=dict(COEFFICIENTS),
        source_parent_context=context, reaction_reference=reference, delta_n_solutes=-1,
        standard_state_conversion_applied_only_after_balanced_source_sum=True,
        missing_microstate_G_is_null=True, ion_pair_populations_inferred=False,
        covalent_quaternization_cycle_modified=False, material_admission=False, source_sha256=sources)
    atomic_write_json(output, result)
    return result
