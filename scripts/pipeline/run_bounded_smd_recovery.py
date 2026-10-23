#!/usr/bin/env python
"""Run the existing source-bound recovery with the separately tested workspace adapter."""
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
import json

ADAPTER = ROOT / "scripts/pipeline/smd_bounded_execution.py"


def compare_previous_pilot(plan, directory):
    previous = Path(plan["source_failed_resource_pilot"]).parent / "energy_and_gradient.json"
    current = directory / "energy_and_gradient.json"
    if not previous.exists():
        return dict(compared=False, reason="Prior OOM prevented a complete fixed-geometry reference")
    a, b = json.loads(previous.read_text()), json.loads(current.read_text())
    differences = dict(maximum_energy_difference_eV=max(abs(a["result"][key]-b["result"][key])
        for key in ("e_smd_eV","e_gas_eV")), gradient_max_difference_Ha_Bohr=abs(a["gradient_max_Ha_Bohr"]-b["gradient_max_Ha_Bohr"]))
    if differences["maximum_energy_difference_eV"] > 1e-5 or differences["gradient_max_difference_Ha_Bohr"] > 1e-6:
        raise ValueError("Bounded fixed-geometry pilot differs from its preserved same-model reference")
    return dict(compared=True, differences=differences, source_sha256={str(p):api.digest(p) for p in (previous,current)})


def checked_controls(plan):
    previous = api.load_module("bounded_control_checker", ROOT / "scripts/pipeline/run_smd_memory_recovery.py")
    previous.ADAPTER = ADAPTER
    evidence = previous.checked_controls(plan)
    for value in plan["SMD_memory_qualification_controls"]:
        path = Path(value).parent / "compact/execution/status.json"
        report = json.loads(path.read_text())
        if (report["status"] != "bounded_execution_finished_requires_control_and_physical_review"
                or report["budget_GiB"] != 80 or report["memory_queries"] <= 0 or report["contracted_R_calls"] <= 0):
            raise ValueError("Bounded execution was not exercised by the control")
        api.verify_sources(report["source_sha256"])
        evidence[str(path)] = api.digest(path)
    return evidence


def main():
    recovery = api.load_module("original_control_gated_SMD_recovery", ROOT / "scripts/pipeline/run_smd_memory_recovery.py")
    recovery.SCRIPT = Path(__file__).resolve()
    recovery.ADAPTER = ADAPTER
    recovery.checked_controls = checked_controls
    def annotate(path, directory, qualification):
        value = json.loads(path.read_text())
        if path.parent.name == "pilot" and path.parent.parent.name == "energy":
            try:
                value["previous_fixed_geometry_comparison"] = compare_previous_pilot(api.read_plan(path.parents[2]),path.parent)
            except ValueError as error:
                value.update(status="failed_preserved",error=str(error))
                atomic_write_json(path,value)
                raise
        value["execution_SMD_memory"] = dict(adapter_sha256=api.digest(recovery.ADAPTER), qualification_controls=qualification,
            workspace_budget_GiB=80, matrix_vector_contraction_reassociated=True,
            scientific_model_changed=False, surface_grid_changed=False, precision_changed=False)
        value["output_sha256"].update({str(p):api.digest(p) for p in directory.rglob("*") if p.is_file()})
        value["output_sha256"].update(qualification)
        atomic_write_json(path,value)
    recovery.annotate = annotate
    recovery.main()


if __name__ == "__main__":
    main()
