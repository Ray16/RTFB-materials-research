from copy import deepcopy
from itertools import islice
import json

import pytest
from rdkit import Chem

from redox.build.attachment_identity import canonical, exact_attachment
from redox.screening import route_stage as routes, streaming_search as search
from redox.workflows import chunked as queue


def case():
    return next(search.family_records())


def task(record):
    result = queue.task_record(record, search.stage_spec("descriptors"))
    result.pop("task_id")
    return result


def evidence(binding):
    return dict(exact_precursor_access_reviewed=True, attachment_selectivity_reviewed=True,
                conditions_defined=True, balanced_route_reviewed=True, redox_core_retention_reviewed=True,
                source_urls=["https://example.org/test-only-route"],
                conditions=dict(solvent="MeCN", reactive_form="test precursor", base="none", counterions="Cl"),
                route_binding=binding)


def screened_payload():
    original = case()
    attachment = exact_attachment(task(original))
    return dict(id="graft_" + attachment["exact_design_key"], fragment=original["payload"]["attachment_fragment"],
                attachment=attachment, descriptor=dict(capacity=999, deltaG_reference=-999), reviewed={}, assessment={})


def test_real_family_recipes_have_balanced_exact_mapped_new_bonds():
    seen = set()
    families = set()
    for record in search.family_records():
        attachment = exact_attachment(task(record))
        binding = attachment["route_binding"]
        molecule = Chem.MolFromSmiles(binding["mapped_product_smiles"])
        indices = {a.GetAtomMapNum(): a.GetIdx() for a in molecule.GetAtoms()}
        bond = molecule.GetBondBetweenAtoms(*[indices[i] for i in binding["new_bond_atom_maps"]])
        assert bond.GetBondType() == Chem.BondType.SINGLE
        assert canonical(binding["mapped_product_smiles"]) == record["identity"]["product_smiles"]
        assert attachment["balanced_reaction"]["atom_balanced"]
        assert attachment["balanced_reaction"]["charge_balanced"]
        seen.add(attachment["exact_design_key"])
        families.add(record["provenance"]["core"])
    assert len(seen) == 837
    assert families == {"viologen", "pmdi", "ndi"}


@pytest.mark.parametrize("field", ["product_smiles", "precursor_smiles", "attachment_fragment"])
def test_changed_attachment_partner_cannot_inherit_old_identity(field):
    record = case()
    record["payload"][field] = "CCO"
    with pytest.raises(ValueError):
        exact_attachment(task(record))


def test_unknown_or_wrong_exact_route_cannot_advance_despite_descriptors():
    payload = screened_payload()
    decision = routes.screen(dict(payload=payload))
    assert decision["route_status"] == "hold"
    assert not decision["proceed_to_full_redox"]
    assert decision["proceed_to_feasibility_model"]
    payload["reviewed"] = evidence(payload["attachment"]["route_binding"])
    assert routes.screen(dict(payload=payload))["proceed_to_full_redox"]
    payload["reviewed"] = deepcopy(payload["reviewed"])
    payload["reviewed"]["route_binding"]["new_bond_atom_maps"][1] += 1
    assert not routes.screen(dict(payload=payload))["proceed_to_full_redox"]


def test_actual_condition_exclusion_and_selection_hold_remain_effective():
    payload = screened_payload()
    payload["reviewed"] = evidence(payload["attachment"]["route_binding"])
    payload["assessment"] = dict(status="on_hold")
    assert not routes.screen(dict(payload=payload))["proceed_to_full_redox"]
    payload["reviewed"]["strongly_unfavorable_under_actual_conditions"] = True
    decision = routes.screen(dict(payload=payload))
    assert decision["route_status"] == "exclude_route"
    assert not decision["proceed_to_feasibility_model"]


def test_ndjson_catalog_spills_to_local_disk_and_rejects_duplicate_ids(tmp_path):
    path = tmp_path / "routes.jsonl"
    row = dict(id="one", evidence=dict(a=1))
    path.write_text(json.dumps(row) + "\n")
    with routes.catalog(path, "evidence") as lookup:
        assert lookup("one") == {"a": 1}
        assert lookup("missing") == {}
    path.write_text((json.dumps(row) + "\n") * 2)
    with pytest.raises(Exception, match="UNIQUE"):
        with routes.catalog(path, "evidence"):
            pass


def test_successor_requires_completed_exact_parent_and_keeps_unknowns(tmp_path):
    parent = tmp_path / "descriptors"
    plan = queue.prepare(islice(search.family_records(), 9), parent, search.stage_spec("descriptors"), chunk_size=20)
    evidence_file, assessment_file = tmp_path / "evidence.json", tmp_path / "assessment.json"
    for path in (evidence_file, assessment_file):
        path.write_text('{"candidates":{}}')
    with pytest.raises(ValueError, match="Complete the descriptor"):
        routes.stage(parent, evidence_file, assessment_file)
    queue.run_chunk(parent, plan["chunks"][0], search.descriptor)
    stage = routes.stage(parent, evidence_file, assessment_file)
    successor = tmp_path / "routes"
    following = queue.prepare(routes.records(parent, evidence_file, assessment_file), successor, stage)
    queue.run_chunk(successor, following["chunks"][0], routes.screen)
    verified = routes.verify(successor)
    assert verified["completed"] == plan["unique_tasks"]
    assert verified["counts"] == {"hold": plan["unique_tasks"]}
    assert verified["quantum_jobs_submitted"] == 0
    with pytest.raises(ValueError, match="exact-route campaign"):
        routes.verify(parent)
