#!/usr/bin/env python3
"""Relocate a prepared review package and reject corrupt/misbound records."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--hessian-python", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source = args.package.resolve()
    base = Path(tempfile.mkdtemp(prefix="redox_review_relocation_"))
    relocated = base / "package"
    shutil.copytree(source, relocated)
    original = str(Path(__file__).resolve().parents[2])
    bootstrap = '''import os,runpy,sys
forbidden=sys.argv[1]
def guard(event,args):
    if event=="open" and isinstance(args[0],(str,bytes)):
        path=os.path.realpath(os.fsdecode(args[0]))
        if path==forbidden or path.startswith(forbidden+os.sep):
            raise RuntimeError("Original checkout/calculation path opened: "+path)
sys.addaudithook(guard)
sys.argv=sys.argv[2:]
runpy.run_path(sys.argv[0],run_name="__main__")
'''
    def run(folder, name, python="/usr/bin/python3", hessian=False):
        target = folder.parent / (name + ".json")
        script = "review_package_hessians.py" if hessian else "reproduce_review_package.py"
        command = [str(python), "-I"] + ([] if hessian else ["-S"])
        command += ["-c", bootstrap, original, str(folder / "code/scripts/release" / script), str(folder), "--out", str(target)]
        env = {"PATH": "/usr/bin:/bin", "LC_ALL": "C.UTF-8"}
        for key in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS"):
            env[key] = "1"
        completed = subprocess.run(command, cwd="/", env=env, capture_output=True, text=True)
        (args.out / (name + ".log")).write_text(completed.stdout + completed.stderr)
        if target.exists():
            shutil.copyfile(target, args.out / (name + ".json"))
        return completed
    positive = run(relocated, "arithmetic")
    if positive.returncode:
        raise ValueError("Relocated offline replay failed; see arithmetic.log")
    negatives = []
    for case in ("changed_bytes", "rehashed_thermal_identity", "rehashed_unbalanced_reaction", "external_path"):
        folder = base / case / "package"
        shutil.copytree(relocated, folder)
        data_path = folder / "data/manifest.json"
        package_path = folder / "package_manifest.json"
        data = json.loads(data_path.read_text())
        package = json.loads(package_path.read_text())
        if case == "changed_bytes":
            file = folder / "data" / data["states"][0]["files"]["energy"]
            file.write_bytes(file.read_bytes() + b" ")
        elif case == "rehashed_thermal_identity":
            relative = next(r["files"]["native"] for r in data["states"] if "native" in r["files"])
            file = folder / "data" / relative
            native = json.loads(file.read_text())
            native["identity"]["geometry_sha256"] = "0" * 64
            file.write_text(json.dumps(native))
            data["files_sha256"][relative] = sha(file)
            package["files_sha256"]["data/" + relative] = sha(file)
        elif case == "rehashed_unbalanced_reaction":
            coefficients = data["reactions"][0]["coefficients"]
            coefficients[next(iter(coefficients))] += 1
        else:
            package["files_sha256"]["../outside-package"] = "0" * 64
        if case in ("rehashed_thermal_identity", "rehashed_unbalanced_reaction"):
            data_path.write_text(json.dumps(data))
            package["files_sha256"]["data/manifest.json"] = sha(data_path)
        package_path.write_text(json.dumps(package))
        failed = run(folder, case)
        if failed.returncode == 0:
            raise ValueError("Negative package test accepted: " + case)
        expected = {"changed_bytes":"artifact differs", "rehashed_thermal_identity":"Thermal identity mismatch",
            "rehashed_unbalanced_reaction":"balance fails", "external_path":"escapes root"}[case]
        if expected not in failed.stderr:
            raise ValueError("Negative test failed for an unrelated reason: " + case)
        negatives.append(dict(case=case, rejected=True, expected_failure=expected))
    if args.hessian_python:
        result = run(relocated, "hessians", args.hessian_python, True)
        if result.returncode:
            raise ValueError("Relocated Hessian reconstruction failed; see hessians.log")
    report = dict(status="relocated_primary_package_checked", checked=time.time(), package=str(source),
        package_manifest_sha256=sha(source / "package_manifest.json"), relocated=str(relocated), cwd="/",
        arithmetic_interpreter="/usr/bin/python3 -I -S", original_checkout_and_calculation_opens_blocked_by_audit_hook=True,
        negative_checks=negatives, optional_Hessian_reconstruction=bool(args.hessian_python),
        source_sha256={str(Path(__file__).resolve()):sha(Path(__file__).resolve())},
        code_and_dataset_publicly_deposited=False, full_GPU_environment_restored=False)
    (args.out / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
