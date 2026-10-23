import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts/pipeline/coordinate_candidate_baseline.py"
spec = importlib.util.spec_from_file_location("candidate_baseline_identity", PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_output_uses_exact_source_key_when_baseline_declaration_omits_it(monkeypatch, tmp_path):
    # Synthetic metadata fixture; no chemical or quantum result is asserted.
    (tmp_path / "plan.json").write_text(json.dumps({"exact_capped_product_key": "EXACT_SOURCE_DESIGN"}))
    (tmp_path / "candidate_snapshot.json").write_text(json.dumps({"id": "fixture"}))
    (tmp_path / "route_snapshot.json").write_text(json.dumps({"gate": {"route_status": "hold", "missing_route_evidence": ["selectivity"]}}))
    (tmp_path / "reaction_terms.json").write_text(json.dumps({"source_sha256": {}, "constant_eV": {}}))
    refs = {str(i): dict(task={"smiles": "C"}, E_smd_eV=i, G_eV=i, state_dir=str(tmp_path)) for i in range(2)}
    monkeypatch.setattr(module.api, "admitted_records", lambda campaign, stage: (refs, []))
    monkeypatch.setattr(module, "load_native_thermal", lambda *args: {"g_thermal_total_eV": 0})
    monkeypatch.setattr(module, "load_record", lambda *args: {"e_smd_eV": 0})
    driver = SimpleNamespace(collect=lambda campaign: {"records": [{"E_eV": None, "G_eV": None}], "checked_energy": 0, "checked_native_G": 0})
    declaration = dict(cases=[dict(campaign=str(tmp_path), occupied_sites=list(sites)) for sites in ((2,), (3,), (2, 3))], shared_empty_reference={"graph_smiles": "C"})
    output = module.collect(tmp_path, declaration, driver)
    assert output["exact_capped_product_key"] == "EXACT_SOURCE_DESIGN"
    assert output["conditional_polymer_cycle"] is None
    assert output["synthesis_admission"] is False


def test_memory_recovery_cannot_change_basis_model_or_task_identity():
    parent = dict(quantum_settings={"method": "fixed"}, thermochemistry={"T": 298.15}, solvent="acetonitrile", tasks=[{"smiles": "exact"}])
    module.same_model(parent, dict(parent))
    for field, bad in (("quantum_settings", {"method": "cheaper"}), ("solvent", "water"), ("tasks", [{"smiles": "another"}])):
        with pytest.raises(ValueError):
            module.same_model(parent, dict(parent, **{field: bad}))
