#!/usr/bin/env python
"""Control-gated SMD temporary-buffer recovery; scientific protocols unchanged."""
import sys
from pathlib import Path
import json

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
ADAPTER = ROOT / "scripts/pipeline/smd_memory_execution.py"


def checked_controls(plan):
    paths = [Path(path) for path in plan["SMD_memory_qualification_controls"]]
    if len(paths) != 2:
        raise ValueError("Both singlet and doublet implementation controls required")
    found = set()
    sources = {}
    for path in paths:
        report = json.loads(path.read_text())
        if report["status"] != "full_original_compact_equivalence_passed" or report["scientific_model_changed"] is not False:
            raise ValueError("SMD memory execution is not qualified")
        found.add(report["mult"])
        api.verify_sources(report["source_sha256"])
        api.verify_sources(report["evidence_sha256"])
        if any(not 0 <= value <= 1e-6 for value in report["differences"].values()):
            raise ValueError("Execution control exceeds its unchanged tolerance")
        if report["source_sha256"].get(str(ADAPTER)) != api.digest(ADAPTER):
            raise ValueError("Adapter differs from the tested implementation")
        sources[str(path)] = api.digest(path)
    if found != {1, 2}:
        raise ValueError("Both spin classes are required")
    return sources


def annotate(path, directory, qualification):
    record = json.loads(path.read_text())
    record["execution_SMD_memory"] = dict(adapter_sha256=api.digest(ADAPTER), qualification_controls=qualification,
        transient_matrix_lifetimes_only=True, scientific_model_changed=False, surface_grid_changed=False)
    record["output_sha256"].update({str(p): api.digest(p) for p in directory.rglob("*") if p.is_file()})
    record["output_sha256"].update(qualification)
    atomic_write_json(path, record)


def main():
    driver = api.load_module("control_gated_memory_recovery", ROOT / "scripts/pipeline/refine_polymer_conformers.py")
    driver.SCRIPT = SCRIPT
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        position = sys.argv.index("--campaign")
        campaign = Path(sys.argv[position + 1]).resolve()
        plan = api.read_plan(campaign)
        qualification = checked_controls(plan)
        memory = api.load_module("qualified_SMD_temporary_buffers", ADAPTER)
        original_energy, original_native = driver.calculate_energy, api.calculate_native_G

        def energy(folder, task, pilot=False):
            directory = Path(plan["raw_root"]) / "execution_SMD_memory" / (task["key"] + ("_resource_pilot" if pilot else "_energy"))
            with memory.compact_smd_execution(directory):
                result = original_energy(folder, task, pilot=pilot)
            path = folder / "energy/pilot/status.json" if pilot else folder / "energy/tasks" / task["key"] / "status.json"
            annotate(path, directory, qualification)
            return result

        def native(folder, record):
            directory = Path(plan["raw_root"]) / "execution_SMD_memory" / (record["task"]["key"] + "_native")
            with memory.compact_smd_execution(directory):
                result = original_native(folder, record)
            annotate(folder / "thermal/tasks" / record["task"]["key"] / "status.json", directory, qualification)
            return result

        driver.calculate_energy = energy
        api.calculate_native_G = native
    driver.main()


if __name__ == "__main__":
    main()
