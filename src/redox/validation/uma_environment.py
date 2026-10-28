"""Inspect both UMA environments without importing Torch or creating a CUDA context."""
import argparse
import json
import os
from pathlib import Path
import subprocess

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm.uma_policy import load_policy
from redox.workflows.environment import worker_environment

PROBE = '''
import importlib.metadata as m, json, pathlib, sys
versions = {}
for name in ("fairchem-core", "torch", "ase", "triton", "nvidia-cutlass-dsl", "cuda-python"):
    try: versions[name] = m.version(name)
    except m.PackageNotFoundError: versions[name] = None
out = dict(python=sys.executable, prefix=sys.prefix, versions=versions)
try:
    d=m.distribution("fairchem-core")
    out["installation_source"] = json.loads(d.read_text("direct_url.json") or "{}")
    p=pathlib.Path(d.locate_file("fairchem/core/models/uma/uma_nitro/_inference.py"))
    out["nitro_executor_exists"] = p.is_file()
    if p.is_file():
        import hashlib
        out["executor_sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
except m.PackageNotFoundError: out["nitro_executor_exists"] = False
print(json.dumps(out))
'''


def inspect_policy(policy, *, verify_model=False):
    report = dict(policy=policy, environments={}, errors=[], cuda_context_created=False,
                  production_qualified=False)
    env = worker_environment(ROOT)
    env["FAIRCHEM_CACHE_DIR"] = policy["fairchem_cache"]
    env["LD_LIBRARY_PATH"] = policy["nitro_library_path"] + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    for name in ("legacy", "nitro"):
        python = policy[name + "_python"]
        try:
            result = subprocess.run([python, "-I", "-c", PROBE], env=env, text=True,
                                    capture_output=True, check=True, timeout=60)
            info = json.loads(result.stdout)
            report["environments"][name] = info
            if not all(info["versions"][key] for key in ("torch", "ase", "fairchem-core")):
                report["errors"].append(f"{name}: missing runtime packages")
            if name == "nitro" and not info["nitro_executor_exists"]:
                report["errors"].append("nitro: standalone executor is missing")
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            report["errors"].append(f"{name}: {exc}")
    if verify_model:
        code = ("import sys,json; sys.path.insert(0," + repr(str(ROOT / "src")) + "); "
                "from redox.qm.uma_model_identity import model_identity; "
                "print(json.dumps(model_identity('uma-s-1p2p1')))")
        try:
            result = subprocess.run([policy["nitro_python"], "-I", "-c", code], env=env,
                                    text=True, capture_output=True, check=True, timeout=180)
            identity = json.loads(result.stdout)
            report["local_model_identity"] = identity
            for key in ("checkpoint_sha256", "reference_sha256"):
                if identity[key] != policy[key]:
                    report["errors"].append(f"Nitro {key} differs from policy")
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            report["errors"].append(f"model: {exc}")
    report["configuration_ok"] = not report["errors"]
    report["status"] = "configuration_verified_numerical_qualification_pending" if not report["errors"] else "configuration_errors"
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-model", action="store_true", help="hash locally cached weights and references; no downloads")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    report = inspect_policy(load_policy(), verify_model=args.verify_model)
    if args.out:
        if args.out.exists():
            parser.error("choose a new output file; existing evidence is preserved")
        atomic_write_json(args.out, report)
    print(json.dumps(report, indent=2))
    return int(not report["configuration_ok"])


if __name__ == "__main__":
    raise SystemExit(main())
