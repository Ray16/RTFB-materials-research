#!/usr/bin/env python
"""Use the standard campaign stages with the controlled XC lifetime adapter."""
from pathlib import Path
from contextlib import contextmanager
import json
import inspect
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def native_host_forecast(surface_points):
    """Use the bounded-PCM driver's formula with a displacement safety margin."""
    if not isinstance(surface_points, int) or surface_points <= 0:
        raise ValueError("Positive measured cavity size required")
    guarded = surface_points + 1024
    return (min(18*guarded*guarded*8, 128*2**30) + 128*2**30) / 2**30


def host_execution_body(body):
    if body.count("else 256") != 2 or body.count("max(budget, required_gib)") != 1:
        raise ValueError("Standard native host-admission template changed")
    body = body.replace("else 256", "else plan.get('native_host_budget_GiB', 256)")
    return body.replace("max(budget, required_gib)",
        "(required_gib if stage == 'thermal' else max(budget, required_gib))")


def main():
    source = ROOT / "scripts/pipeline/workspace_polymer_campaign.py"
    module = api.load_module("XC_workspace_campaign_stages", source)
    module.SCRIPT = Path(__file__).resolve()
    module.ADAPTER = ROOT / "scripts/pipeline/xc_workspace_execution.py"
    module.CONTROL_CHECKER = ROOT / "scripts/validation/qualify_xc_checkpoint.py"
    checker = api.load_module("XC_workspace_campaign_controls", module.CONTROL_CHECKER)

    def controls(plan):
        paths = [Path(p) for p in plan["SMD_memory_qualification_controls"]]
        roots = {p.parents[2] for p in paths}
        if len(paths) != 2 or len(roots) != 1:
            raise ValueError("Matched XC execution controls required")
        return checker.checked_controls(roots.pop(), plan["execution_workspace_GiB"])

    module.controls = controls
    body = host_execution_body(inspect.getsource(module.execution))
    namespace = dict(vars(module))
    exec(compile(body, str(Path(__file__).resolve())+":host_admission", "exec"), namespace)
    module.execution = namespace["execution"]
    if len(sys.argv) > 1 and sys.argv[1] == "worker":
        original_loader = api.load_module

        def loader(name, path):
            loaded = original_loader(name, path)
            if name != "workspace_campaign_adapter":
                return loaded
            original_context = loaded.compact_smd_execution

            @contextmanager
            def observed(directory, budget_GiB):
                with original_context(directory, budget_GiB=budget_GiB) as proof:
                    import cupy
                    from gpu4pyscf.grad import rks
                    namespace = rks.get_exc_full_response.__globals__
                    available = namespace["get_avail_mem"]

                    def logged(exclude_memory_pool=False):
                        value = available(exclude_memory_pool)
                        free, total = cupy.cuda.runtime.memGetInfo()
                        pool = cupy.get_default_memory_pool()
                        row = dict(time=time.time(), returned_available_bytes=value,
                            GPU_total_bytes=total, GPU_free_bytes=free,
                            pool_total_bytes=pool.total_bytes(), pool_used_bytes=pool.used_bytes())
                        if value <= 0:
                            try:
                                row["resident_GPU_processes"] = subprocess.check_output([
                                    "nvidia-smi", "--query-compute-apps=pid,gpu_uuid,used_gpu_memory",
                                    "--format=csv,noheader"], text=True, timeout=10)
                            except (OSError, subprocess.SubprocessError) as error:
                                row["observation_error"] = repr(error)
                        with (Path(directory) / "RKS_available_memory.jsonl").open("a") as stream:
                            stream.write(json.dumps(row)+"\n")
                        return value

                    namespace["get_avail_mem"] = logged
                    try:
                        yield proof
                    finally:
                        namespace["get_avail_mem"] = available

            loaded.compact_smd_execution = observed
            return loaded

        api.load_module = loader
    module.main()
    if sys.argv[1] == "launch" and "--stage" in sys.argv and sys.argv[sys.argv.index("--stage")+1] == "thermal":
        campaign = Path(sys.argv[sys.argv.index("--campaign")+1]).resolve()
        plan = api.read_plan(campaign)
        import socket
        p = campaign / "thermal" / (socket.gethostname().split('.')[0]+"_RAM_admission.json")
        if p.exists() and "native_host_budget_GiB" in plan:
            data = json.loads(p.read_text())
            data["legacy_default_host_budget_GiB"] = data["worker_budget_GiB"]
            data["worker_budget_GiB"] = plan["native_host_budget_GiB"]
            data["forecast_scope"] = "Measured cavity plus displacement margin; original curvature/cache/headroom formula"
            atomic_write_json(p, data)
    if sys.argv[1].startswith("prepare-"):
        out = Path(sys.argv[sys.argv.index("--out") + 1]).resolve()
        plan = api.read_plan(out)
        plan["producer_sha256"][str(source)] = api.digest(source)
        plan["XC_temporary_lifetime_adapter"] = str(module.ADAPTER)
        snapshot = out / "host_admission_execution.py"
        snapshot.write_text(body)
        plan["producer_sha256"][str(snapshot)] = api.digest(snapshot)
        if sys.argv[1] == "prepare-native":
            parent = Path(sys.argv[sys.argv.index("--source") + 1]).resolve()
            report, _ = module.checked_energy(parent)
            evidence = dict(report["output_sha256"])
            workflow = parent / "workspace_energy_status.json"
            if workflow.exists():
                evidence.update(json.loads(workflow.read_text()).get("execution_sha256", {}))
            build_paths = [Path(p) for p in evidence if p.endswith("build_execution/status.json")]
            sizes = []
            for p in build_paths:
                if api.digest(p) != evidence[str(p)]:
                    raise ValueError("Cavity-size source evidence changed")
                sizes += [r["surface_points"] for r in json.loads(p.read_text())["builds"]
                          if r["atoms"] == report["task"]["atoms"]]
                plan["input_sha256"][str(p)] = api.digest(p)
            if sizes:
                plan["native_host_budget_GiB"] = native_host_forecast(max(sizes))
                plan["host_required_bytes"] = int(plan["native_host_budget_GiB"]*2**30)
                plan["native_host_forecast"] = dict(measured_max_surface_points=max(sizes),
                    displacement_margin_points=1024, unchanged_cache_limit_GiB=128,
                    unchanged_host_headroom_GiB=128, source_build_reports=list(map(str, build_paths)))
        atomic_write_json(out / "plan.json", plan)


if __name__ == "__main__":
    main()
