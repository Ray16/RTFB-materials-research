import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]


def script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts/pipeline" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


prepare = script("prepare_crest_refinement")
refine = script("refine_polymer_conformers")
coordinate = script("coordinate_crest_refinement")


def test_plain_xyz_units_are_not_inferred_from_an_arbitrary_number():
    assert prepare.plain_xyz_energy("1\n -2.1234\nH 0 0 0\n") == -2.1234
    for comment in ["energy=-2.1234", "Etot= -2.1234 g norm=0.1", "energy=1 energy_units=eV", "nan", "inf"]:
        with pytest.raises(ValueError):
            prepare.plain_xyz_energy("1\n" + comment + "\nH 0 0 0\n")


def test_running_search_never_admits_its_intermediate_frames(tmp_path):
    raw = tmp_path / "search"
    raw.mkdir()
    (raw / "status.json").write_text(json.dumps({"status": "running"}))
    (raw / "crest_conformers.xyz").write_text("1\n-1\nH 0 0 0\n")
    assert prepare.source_records({"source_search": str(raw)}) is None


def test_final_label_without_terminal_evidence_is_rejected(tmp_path):
    (tmp_path / "status.json").write_text(json.dumps({"status": "low_level_ensemble_complete_DFT_pending"}))
    with pytest.raises(ValueError, match="terminal"):
        prepare.source_records({"source_search": str(tmp_path)})


def test_two_configurations_are_kept_distinct():
    assert prepare.configuration({"backbone_stereo_tags": ["CW"] * 4}) == "all_same"
    assert prepare.configuration({"backbone_stereo_tags": ["CW", "CCW", "CW", "CCW"]}) == "alternating"
    with pytest.raises(ValueError):
        prepare.configuration({"backbone_stereo_tags": ["CW"]})


def rows(n):
    return [dict(key=f"c{i:02d}", E_eV=-100.0, G_eV=-99.0) for i in range(n)]


def test_nested_batches_do_not_count_duplicate_discoveries_as_entropy():
    data = refine.finite_batches(rows(20), np.zeros((20, 20)))
    for value in data.values():
        assert value["10"]["unique_minima"] == value["20"]["unique_minima"] == 1
        assert value["convergence"]["ensemble_G_change_kJ_mol"] == 0
        assert value["convergence"]["global_ensemble_convergence_proven"] is False


def test_later_duplicate_representative_cannot_leak_into_earlier_batch():
    data = rows(20)
    data[19]["G_eV"] = -101.0
    result = refine.finite_batches(data, np.zeros((20, 20)))["0.1"]
    assert result["10"]["G_eV"] == -99
    assert result["20"]["G_eV"] == -101
    assert result["convergence"]["within_batch_tolerance"] is False


def test_nineteen_starts_cannot_satisfy_the_twenty_start_comparison():
    result = refine.finite_batches(rows(19), np.zeros((19, 19)))["0.1"]
    assert result["20"]["G_eV"] is None
    assert result["convergence"]["within_batch_tolerance"] is None


def test_missing_native_G_stays_missing_at_all_duplicate_cutoffs():
    data = rows(20)
    data[15]["G_eV"] = None
    result = refine.finite_batches(data, np.zeros((20, 20)))
    assert all(value["20"]["G_eV"] is None for value in result.values())


def test_native_dispatch_waits_for_every_declared_energy_parent(monkeypatch, tmp_path):
    monkeypatch.setattr(refine.api, "admitted_records", lambda folder, stage: ({}, []))
    assert refine.backend(tmp_path, "thermal").phase_tasks(tmp_path, {"tasks": [{}, {}]}, "thermal") == []


def test_hessian_startup_promises_prevent_memory_overcommit():
    gib = 1024**3
    promises = [dict(budget_bytes=256 * gib, observed_unique_bytes=1 * gib)] * 3
    assert refine.additional_host_slots(1024 * gib, promises, 256) == 0
    assert refine.additional_host_slots(2048 * gib, promises, 256) == 3


def test_already_allocated_unique_RAM_is_not_reserved_twice():
    gib = 1024**3
    promises = [dict(budget_bytes=256 * gib, observed_unique_bytes=200 * gib)] * 2
    assert refine.additional_host_slots(1000 * gib, promises, 256) == 2


def test_staging_maps_core_indices_into_reordered_XYZ(monkeypatch, tmp_path):
    # Synthetic producer fixture only: an embedded one-repeat graph exercises
    # staging/atom mapping. No CREST completion or quantum result is asserted.
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase.io import read, write
    from redox.build.polymer_cases import build_case
    from redox.build.structures import to_xyz
    from redox.build.polymer_contacts import pi_core_groups

    case = build_case(prepare.FRAGMENT, 1, [1], [1], [])
    case.update(key="fixture", mult=1)
    mol = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
    assert AllChem.EmbedMolecule(mol, randomSeed=47) == 0
    source = tmp_path / "embedded.xyz"
    source.write_text(to_xyz(mol, "synthetic geometry fixture"))
    atoms = read(source)
    # Put hydrogens before heavy atoms. Treating graph indices as XYZ indices
    # would now select the wrong atoms for the conjugated core.
    write(source, atoms[::-1], format="xyz")
    row = dict(index=0, path=str(source), sha256=prepare.api.digest(source), GFN2_energy_Ha=-100.0)
    monkeypatch.setattr(prepare, "source_records", lambda entry: (case, [row], [], {}))
    root = tmp_path / "cohort"
    root.mkdir()
    request = dict(maximum_initial_starts=20, ordered_DFT_batches=[10, 20], raw_root=str(tmp_path / "raw"))
    (root / "request.json").write_text(json.dumps(request))
    entry = dict(key="fixture", fragment=prepare.FRAGMENT, motif_id="fixture_only", role="method_control",
        configuration="fixture_only", source_search="synthetic_fixture", mandatory_references=[])
    result = prepare.stage_case(root, request, entry)
    assert result["starts"] == 1
    task = prepare.api.read_plan(root / "fixture")["tasks"][0]
    geometry = refine.api.load_module("fixture_mapping", ROOT / "scripts/pipeline/polymer_conformer_thermals.py")
    # Rebuild the declared atom order from XYZ; the production core set must
    # contain carbon/nitrogen ring atoms, never the leading XYZ hydrogens.
    xyz = root / "fixture/live_inputs" / task["key"] / "input.xyz"
    observed = read(xyz)
    selected = task["pi_groups"][0]
    assert len(selected) == len(pi_core_groups(mol, prepare.FRAGMENT)[0])
    assert all(observed[i].symbol in ("C", "N") for i in selected)
    assert geometry.alignment_molecule(case["smiles"], xyz).GetNumAtoms() == Chem.RemoveHs(mol).GetNumAtoms()


def cycle_fixture(tmp_path):
    cases = []
    reports = {}
    values = {(): 0, (2,): 1, (3,): 2, (4,): 3, (2, 3): 4, (2, 4): 5}
    for configuration, factor in (("all_same", 1), ("alternating", 10)):
        for i, (occupied, value) in enumerate(values.items()):
            key = configuration + str(i)
            case = dict(occupied_sites=list(occupied), smiles="fixture_graph_" + key,
                repeat_units=5, reactive_sites=[2, 3, 4], backbone_stereo_fingerprint=[configuration])
            path = tmp_path / (key + ".json")
            path.write_text(json.dumps(case))
            cases.append(dict(configuration=configuration, key=key, case_path=str(path)))
            reports[key] = dict(primary_G_eV=-1000 + value * factor, partition_protocol={"fixture_protocol": True},
                temperature_K=298.15, solvent="acetonitrile", case=case)
    return dict(cases=cases, temperature_K=298.15, solvent="acetonitrile"), reports


def test_conditional_cycles_keep_backbone_configurations_separate(tmp_path):
    from redox.properties.conformer_ensemble import EV_KJ_MOL
    request, reports = cycle_fixture(tmp_path)
    cycles = coordinate.matched_cycles(request, reports)
    assert len(cycles) == 4
    for cycle in cycles:
        factor = 1 if cycle["configuration"] == "all_same" else 10
        assert cycle["deltaDeltaG_kJ_mol"] == pytest.approx(factor * EV_KJ_MOL)


@pytest.mark.parametrize("corruption", ["protocol", "temperature", "identity"])
def test_cycle_join_rejects_missing_context_or_wrong_graph(tmp_path, corruption):
    request, reports = cycle_fixture(tmp_path)
    row = next(iter(reports.values()))
    if corruption == "protocol":
        row["partition_protocol"] = None
    elif corruption == "temperature":
        row["temperature_K"] = 310
    else:
        row["case"] = dict(row["case"], smiles="another_graph")
    with pytest.raises(ValueError):
        coordinate.matched_cycles(request, reports)


def test_live_launch_observer_is_repolled_without_duplicate_submission(tmp_path):
    # A short stdlib-only child models a launcher, not a quantum calculation.
    handle = tmp_path / "handle.json"
    outcome = tmp_path / "outcome.json"
    child = "import os,json,time,sys; from pathlib import Path; time.sleep(.5); Path(sys.argv[1]).write_text(json.dumps(dict(pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],status='launch_observation_completed_outputs_preserved',finished=time.time())))"
    payload = dict(root=str(tmp_path), handle=str(handle), outcome=str(outcome), log=str(tmp_path / "launch.log"), argv=[sys.executable, "-c", child, str(outcome)])
    def invoke():
        result = subprocess.run([sys.executable, "-c", coordinate.LAUNCH_REMOTE], input=json.dumps(payload), capture_output=True, text=True, check=True, timeout=5)
        return json.loads(result.stdout)
    first = invoke()
    second = invoke()
    assert first["status"] == "launch_observer_started_scientific_starts_unconfirmed"
    assert second["status"] == "same_launch_observer_live"
    assert first["handle"]["pid"] == second["handle"]["pid"]
    deadline = time.time() + 3
    while not outcome.exists() and time.time() < deadline:
        time.sleep(.05)
    assert outcome.exists()
    # Losing the observation outcome is not permission to launch a duplicate.
    time.sleep(.1)
    outcome.unlink()
    third = invoke()
    assert third["status"] == "ended_launch_outcome_unconfirmed_preserve"
    assert third["handle"]["pid"] == first["handle"]["pid"]


def test_unverified_device_is_skipped_without_a_task_failure(monkeypatch, tmp_path):
    from redox.workflows import resources
    def unavailable(gpu):
        raise subprocess.TimeoutExpired(["nvidia-smi"], 15)
    monkeypatch.setattr(resources, "gpu_suitability", unavailable)
    class FakeBackend:
        def launch(self, campaign, stage, pilot):
            with pytest.raises(RuntimeError, match="no task attempt"):
                resources.gpu_suitability(0)
    monkeypatch.setattr(refine, "backend", lambda campaign, stage: FakeBackend())
    observations = refine.launch(tmp_path, "energy", False)
    assert observations[0]["status"] == "identity_temporarily_unverified"
    assert resources.gpu_suitability is unavailable
    assert not list(tmp_path.glob("**/failures/*.json"))
