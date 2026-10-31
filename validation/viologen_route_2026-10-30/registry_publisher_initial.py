"""Consolidate exact local-graft designs while preserving every source declaration."""
from collections import Counter
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from redox.build.structures import decorate
from redox.build.grafting import precursor
from redox.core.io import atomic_write_json


def merge_sources(sources):
    records, aliases = {}, {}
    for source, role, rows in sources:
        for index, original in enumerate(rows):
            row = original.get('group', original)
            fragment = row['frag']
            mol = decorate('Cc1ccc(C[*:1])cc1', fragment)
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            identifier = row['id']
            declared = (original.get('grafted_smiles') or row.get('exact_grafted_smiles')
                        or row.get('product_smiles') or row.get('route_binding', {}).get('product_smiles'))
            if declared and Chem.MolToSmiles(Chem.MolFromSmiles(declared), isomericSmiles=True) != smiles:
                raise ValueError(f'{identifier}: declared product differs from attachment graph')
            key = hashlib.sha256(smiles.encode()).hexdigest()
            if identifier in aliases and aliases[identifier] != key:
                raise ValueError(f'Alias {identifier} belongs to different chemical graphs')
            aliases[identifier] = key
            pre, kind = precursor(fragment)
            states = row.get('states', row.get('proposed_charge_spin_states'))
            if isinstance(states, str):
                states = json.loads(states)
            record = records.setdefault(key, dict(graph_key=key, product_smiles=smiles,
                inchikey=Chem.MolToInchiKey(mol), aliases=[], families=[], roles=[], declarations=[],
                evaluations=[], recipes=[], route_status='hold', material_decision='pending'))
            for field, value in [('aliases', identifier), ('families', row.get('family', original.get('core', 'unclassified'))), ('roles', role)]:
                if value not in record[field]: record[field].append(value)
            record['declarations'].append(dict(source=source, row=index, id=identifier,
                fragment=fragment, precursor_smiles=Chem.MolToSmiles(pre), reaction_class=kind,
                states=states, state_scope='declared hypothesis; field-specific QC required',
                design_factors=row.get('design_factors', {})))
            for recipe in original.get('recipes', []):
                if recipe not in record['recipes']: record['recipes'].append(recipe)
            if role in ('evaluated_candidate', 'comparison_control'):
                if row.get('parent_metrics_inherited') is True:
                    raise ValueError('Parent properties cannot enter exact-graft evaluations')
                record['evaluations'].append(dict(source=source, row=index,
                    reported_values={k: v for k, v in row.items()
                        if k in ('E1_vs_Fc_V', 'E2_vs_Fc_V', 'wave_spacing_V',
                                 'dG_disproportionation_kJ_mol', 'lambda_i_eV',
                                 'lambda_outer_eV', 'capacity_nominal_mAh_g',
                                 'accessible_capacity_mAh_g', 'MW')},
                    coverage={k: v for k, v in row.items() if k.endswith('_status') or k.endswith('_QC_states')},
                    raw_records_revalidated=False))
    return sorted(records.values(), key=lambda r: r['graph_key'])


def publish(spec_path, archive, output, *, route_evidence=None):
    spec_path, archive, output = Path(spec_path), Path(archive).resolve(), Path(output)
    if output.exists(): raise ValueError('Use a fresh publication directory')
    spec = json.loads(spec_path.read_text())
    sources, inventory, frozen = [], [], {}
    for item in spec['sources']:
        base = spec_path.resolve().parents[1] if item.get('root') == 'checkout' else archive
        path = base / item['path']
        raw = path.read_bytes(); digest = hashlib.sha256(raw).hexdigest()
        payload = json.loads(raw)
        rows = payload[item['section']] if item.get('section') else payload
        name = item['name']
        if name in frozen: raise ValueError('Duplicate source name')
        frozen[name] = raw
        inventory.append(dict(name=name, original_path=str(path), sha256=digest,
                              snapshot='inputs/'+name+'.json', section=item.get('section'), role=item['role'], rows=len(rows)))
        sources.append((name, item['role'], rows))
    records = merge_sources(sources)
    # Route decisions are computed against each actual alias/fragment; a pass
    # on one alias is retained with its evidence rather than inferred for others.
    import importlib.util
    gate_path = Path(__file__).resolve().parents[3] / 'scripts/pipeline/screen_grafting_routes.py'
    module = importlib.util.spec_from_file_location('unified_route_gate', gate_path)
    gate = importlib.util.module_from_spec(module); module.loader.exec_module(gate)
    evidence_path = Path(route_evidence).resolve() if route_evidence is not None else archive / spec['route_evidence']
    assessment_path = archive / spec['assessment']
    route_bytes = evidence_path.read_bytes()
    assessment_bytes = assessment_path.read_bytes()
    evidence = json.loads(route_bytes)['candidates']
    assessment = json.loads(assessment_bytes)
    for record in records:
        decisions = [gate.screen(dict(id=d['id'], frag=d['fragment'], product_smiles=record['product_smiles']), evidence, assessment)
                     for d in record['declarations']]
        record['route_decisions'] = [dict(id=d['id'], status=d['route_status'],
            missing=d.get('missing_route_evidence', []), reasons=d['reasons']) for d in decisions]
        statuses = {d['route_status'] for d in decisions}
        # Conflicting route-specific exclusions require review; no silent override.
        record['route_status'] = ('review_conflict' if 'proceed' in statuses and any(s.startswith('exclude') for s in statuses)
                                  else 'proceed' if 'proceed' in statuses else 'hold')
    for item in inventory:
        if hashlib.sha256(Path(item['original_path']).read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('Input changed during registry assembly')
    if evidence_path.read_bytes() != route_bytes or assessment_path.read_bytes() != assessment_bytes:
        raise ValueError('Route evidence or assessment changed during registry assembly')
    output.mkdir(parents=True)
    (output/'inputs').mkdir()
    for name, raw in frozen.items(): (output/'inputs'/f'{name}.json').write_bytes(raw)
    for name, path in [('route_evidence', evidence_path), ('assessment', assessment_path), ('spec', spec_path)]:
        raw=path.read_bytes(); (output/'inputs'/f'{name}.json').write_bytes(raw)
        inventory.append(dict(name=name, original_path=str(path.resolve()), sha256=hashlib.sha256(raw).hexdigest(), snapshot='inputs/'+name+'.json'))
    atomic_write_json(output/'registry.json', dict(schema=1, candidates=records))
    summary=dict(exact_graphs=len(records), aliases=len({a for r in records for a in r['aliases']}),
        declarations=sum(len(r['declarations']) for r in records),
        evaluated_graphs=sum(bool(r['evaluations']) for r in records),
        route_counts=dict(Counter(r['route_status'] for r in records)),
        full_material_admissions=0, raw_quantum_records_revalidated=False,
        new_quantum_jobs=0, inputs=inventory,
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        route_gate_sha256=hashlib.sha256(gate_path.read_bytes()).hexdigest(),
        property_validation_scope='Exact-graph joins and recorded field coverage; raw QC and external accuracy require separate evidence')
    atomic_write_json(output/'summary.json', summary)
    return summary
