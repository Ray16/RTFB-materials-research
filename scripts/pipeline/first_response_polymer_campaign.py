#!/usr/bin/env python
"""Native polymer stages with independently checked blocked PCM response.

Preparation requires the exact geometry's response/resource qualification.
Full native resource and physical gates still run before any G is accepted.
"""
import inspect
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/smd_first_response_workspace_v2.py"
CHECKER = ROOT / "scripts/validation/qualify_pcm_first_response_v2.py"


def response_evidence(folder, report):
    status_path, plan_path, review_path = (folder / n for n in
        ("status.json", "plan.json", "independent_review.json"))
    status, plan, review = (json.loads(p.read_text()) for p in (status_path, plan_path, review_path))
    if (status["status"] != "response_resource_and_sampled_derivatives_passed"
            or review["status"] != "fixed_density_response_matrices_independently_reconstructed"
            or review["maximum_relative_error"] > .01
            or review["qualified_native_Hessian"] or review["qualified_G"]
            or review["observed_peak_GPU_MiB"] > 32768):
        raise ValueError("The isolated response/independent review is not qualified within scope")
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    api.verify_sources(status["output_sha256"])
    api.verify_sources(review["source_sha256"])
    if (api.digest(plan["input_xyz"]) != api.digest(Path(report["state_dir"]) / "opt.xyz")
            or plan["task"]["key"] != report["task"]["key"]
            or review["source_sha256"].get(str(status_path)) != api.digest(status_path)):
        raise ValueError("The response qualification belongs to another geometry or state")
    return {str(p): api.digest(p) for p in (status_path, plan_path, review_path)}


def main():
    standard_path = ROOT / "scripts/pipeline/workspace_polymer_campaign.py"
    host_path = ROOT / "scripts/pipeline/xc_workspace_campaign.py"
    alias_path = ROOT / "scripts/pipeline/materialize_imported_energy_parent.py"
    module = api.load_module("first_response_standard_campaign", standard_path)
    host = api.load_module("first_response_existing_host_forecast", host_path)
    module.SCRIPT = Path(__file__).resolve()
    module.ADAPTER, module.CONTROL_CHECKER = ADAPTER, CHECKER
    body = host.host_execution_body(inspect.getsource(module.execution))
    namespace = dict(vars(module))
    exec(compile(body, str(Path(__file__).resolve())+":host_admission", "exec"), namespace)
    module.execution = namespace["execution"]
    witnesses = None
    if sys.argv[1] == "prepare-native":
        if "--response-proof" not in sys.argv:
            raise ValueError("Native preparation requires the exact-geometry response proof")
        index = sys.argv.index("--response-proof")
        folder = Path(sys.argv[index+1]).resolve()
        del sys.argv[index:index+2]
        source = Path(sys.argv[sys.argv.index("--source")+1]).resolve()
        report, _ = module.checked_energy(source)
        witnesses = response_evidence(folder, report)
    module.main()
    if sys.argv[1] == "prepare-native":
        out = Path(sys.argv[sys.argv.index("--out")+1]).resolve()
        plan = api.read_plan(out)
        plan["input_sha256"].update(witnesses)
        plan["response_qualification"] = str(folder)
        sources = (Path(__file__).resolve(), standard_path, host_path, alias_path)
        plan["producer_sha256"].update({str(p): api.digest(p) for p in sources})
        evidence = dict(report["output_sha256"])
        direct = source / "workspace_energy_status.json"
        if direct.exists():
            evidence.update(json.loads(direct.read_text()).get("execution_sha256", {}))
        sizes = []
        for path, expected in evidence.items():
            if not path.endswith("build_execution/status.json"):
                continue
            if api.digest(path) != expected:
                raise ValueError("Measured cavity backing changed")
            sizes.extend(row["surface_points"] for row in json.loads(Path(path).read_text())["builds"]
                         if row["atoms"] == report["task"]["atoms"])
            plan["input_sha256"][path] = expected
        if not sizes:
            raise ValueError("Measured parent cavity size is required")
        plan["native_host_budget_GiB"] = host.native_host_forecast(max(sizes))
        plan["host_required_bytes"] = int(plan["native_host_budget_GiB"]*2**30)
        plan["native_host_forecast"] = dict(measured_max_surface_points=max(sizes),
            displacement_margin_points=1024, unchanged_cache_limit_GiB=128,
            unchanged_host_headroom_GiB=128)
        snapshot = out / "host_admission_execution.py"
        snapshot.write_text(body)
        plan["producer_sha256"][str(snapshot)] = api.digest(snapshot)
        aliases = api.load_module("first_response_parent_aliases", alias_path).materialize(out)
        for row in aliases["aliases"]:
            plan["input_sha256"][row["alias"]] = row["sha256"]
        atomic_write_json(out / "plan.json", plan)


if __name__ == "__main__":
    main()
