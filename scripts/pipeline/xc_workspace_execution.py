"""Release the completed density-loop AO batch before XC grid-response work.

This changes only temporary lifetimes and synchronization. The rho, XC/grid
equations, integration points and existing workspace budget are retained.
"""
from contextlib import contextmanager
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import sys
import subprocess
import time


def release_density_batch(source):
    old = "    assert g1 == ngrids\n\n    exc, vxc ="
    new = ("    assert g1 == ngrids\n"
           "    _workspace_trace('density_loop_complete')\n"
           "    cupy.cuda.get_current_stream().synchronize()\n"
           "    del ao, idx, weight, dms_masked\n"
           "    _workspace_trace('density_batch_released')\n\n"
           "    exc, vxc =")
    if source.count(old) != 1:
        raise ValueError("Installed RKS full-response density-loop template changed")
    return source.replace(old, new)


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    import cupy
    from gpu4pyscf.grad import rks
    source = Path(__file__).with_name("smd_gradient_workspace.py")
    spec = importlib.util.spec_from_file_location("qualified_gradient_workspace_underlay", source)
    previous = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(previous)
    directory = Path(directory)
    with previous.compact_smd_execution(directory, budget_GiB=budget_GiB) as proof:
        original = rks.get_exc_full_response
        before = inspect.getsource(original)
        after = release_density_batch(before)
        original_path = directory / "original_RKS_XC_response.py"
        replacement_path = directory / "released_RKS_XC_response.py"
        header = "\n".join(Path(rks.__file__).read_text().splitlines()[:13])+"\n\n"
        original_path.write_text(header+before)
        replacement_path.write_text(header+after)
        proof.update(XC_density_batch_releases=0, XC_equations_changed=False)
        sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        proof["source_sha256"].update({str(p): sha(p) for p in
            (Path(__file__).resolve(), source, Path(rks.__file__), original_path, replacement_path)})

        def trace(stage):
            free, total = cupy.cuda.runtime.memGetInfo()
            pool = cupy.get_default_memory_pool()
            row = dict(stage=stage, time=time.time(), GPU_total_bytes=total,
                GPU_free_bytes=free, pool_total_bytes=pool.total_bytes(),
                pool_used_bytes=pool.used_bytes(), pool_free_bytes=pool.free_bytes(),
                available_bytes=rks.get_avail_mem())
            if row["available_bytes"] <= 0:
                try:
                    row["resident_GPU_processes"] = subprocess.check_output([
                        "nvidia-smi", "--query-compute-apps=pid,gpu_uuid,used_gpu_memory",
                        "--format=csv,noheader"], text=True, timeout=10)
                except (OSError, subprocess.SubprocessError) as error:
                    row["resident_process_observation_error"] = repr(error)
            if stage == "density_batch_released":
                proof["XC_density_batch_releases"] += 1
            with (directory / "XC_workspace_trace.jsonl").open("a") as stream:
                stream.write(json.dumps(row)+"\n")

        namespace = dict(original.__globals__, _workspace_trace=trace)
        local = {}
        exec(compile(after, str(replacement_path), "exec"), namespace, local)
        replacement = local["get_exc_full_response"]
        for name, module in list(sys.modules.items()):
            if name.startswith("gpu4pyscf") and module is not None:
                for key, value in list(vars(module).items()):
                    if value is original:
                        setattr(module, key, replacement)
        try:
            yield proof
        finally:
            for name, module in list(sys.modules.items()):
                if name.startswith("gpu4pyscf") and module is not None:
                    for key, value in list(vars(module).items()):
                        if value is replacement:
                            setattr(module, key, original)
