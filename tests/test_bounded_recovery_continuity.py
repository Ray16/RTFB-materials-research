"""Same-model memory recovery must preserve available larger-case references."""
import importlib.util
import json
from pathlib import Path
import pytest

SPEC=importlib.util.spec_from_file_location("bounded_recovery",Path(__file__).resolve().parents[1]/"scripts/pipeline/run_bounded_smd_recovery.py")
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_no_reference_is_reported_as_uncompared(tmp_path):
    assert MODULE.compare_previous_pilot({"source_failed_resource_pilot":str(tmp_path/"missing/status.json")},tmp_path)["compared"] is False


def test_changed_energy_blocks_recovery(tmp_path):
    old,new=tmp_path/"old",tmp_path/"new";old.mkdir();new.mkdir()
    record=dict(result={"e_smd_eV":-100.,"e_gas_eV":-99.},gradient_max_Ha_Bohr=.01)
    (old/"energy_and_gradient.json").write_text(json.dumps(record));(new/"energy_and_gradient.json").write_text(json.dumps(record))
    plan={"source_failed_resource_pilot":str(old/"status.json")}
    assert MODULE.compare_previous_pilot(plan,new)["compared"] is True
    record["result"]["e_smd_eV"]+=.1;(new/"energy_and_gradient.json").write_text(json.dumps(record))
    with pytest.raises(ValueError,match="differs"):
        MODULE.compare_previous_pilot(plan,new)
