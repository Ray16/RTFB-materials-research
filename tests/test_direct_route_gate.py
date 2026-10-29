"""The standalone gate must bind route reviews to the exact attachment."""
from copy import deepcopy
import importlib.util
from pathlib import Path

from redox.build.attachment_identity import exact_attachment
from redox.screening import streaming_search as search
from redox.workflows import chunked as queue

spec = importlib.util.spec_from_file_location('direct_route_gate', Path(__file__).resolve().parents[1] / 'scripts/pipeline/screen_grafting_routes.py')
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def test_direct_gate_requires_matching_mapped_bond():
    record = next(search.family_records())
    task = queue.task_record(record, search.stage_spec('descriptors'))
    task.pop('task_id')
    binding = exact_attachment(task)['route_binding']
    row = {'id': 'reviewed_id', 'frag': record['payload']['attachment_fragment']}
    review = dict(exact_precursor_access_reviewed=True, attachment_selectivity_reviewed=True,
                  conditions_defined=True, balanced_route_reviewed=True,
                  redox_core_retention_reviewed=True, source_urls=['https://example.org/test'],
                  conditions=dict(solvent='MeCN', reactive_form='test', base='none', counterions='Cl'),
                  route_binding=binding)
    def decision():
        return gate.screen(row, {row['id']: review}, {'candidates': {}})
    assert decision()['proceed_to_full_redox']
    review['route_binding'] = deepcopy(binding)
    review['route_binding']['new_bond_atom_maps'][1] += 1
    assert not decision()['proceed_to_full_redox']
    review.pop('route_binding')
    assert not decision()['proceed_to_full_redox']
    review['route_binding'] = binding
    row['product_smiles'] = 'CCO'
    assert not decision()['proceed_to_full_redox']
