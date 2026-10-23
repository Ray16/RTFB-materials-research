#!/usr/bin/env python
"""Isolated fixed-geometry footprint qualification; never a candidate G result.

Inventory both screened DF tensor shapes before allocating their host storage.
Run the complete existing energy/gradient pilot only when host RAM is sufficient.
"""
import argparse
import gc
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/smd_workspace_execution.py"


def checked_controls(paths,budget):
    previous=api.load_module("workspace_control_checker",ROOT/"scripts/pipeline/run_smd_memory_recovery.py")
    previous.ADAPTER=ADAPTER
    evidence=previous.checked_controls({"SMD_memory_qualification_controls":paths})
    for value in paths:
        directory=Path(value).parent/"compact/execution"
        for path in (directory/"status.json",directory/"bounded/status.json"):
            data=json.loads(path.read_text())
            if data["budget_GiB"]!=budget or "finished_requires_control_and_physical_review" not in data["status"]:
                raise ValueError("Workspace control budget/status differs")
            api.verify_sources(data["source_sha256"]);evidence[str(path)]=api.digest(path)
        profile=json.loads((directory/"status.json").read_text())
        if profile["prebuilt_solvent_calls"]<=0 or not profile["SCF_calls"]:
            raise ValueError("Allocation order/checkpoint path not exercised")
        bounded=json.loads((directory/"bounded/status.json").read_text())
        if bounded["memory_queries"]<=0 or bounded["contracted_R_calls"]<=0:
            raise ValueError("Bounded workspace/response contraction not exercised")
    return evidence


def prepare(args):
    parent=args.parent.resolve();out=args.campaign.resolve();old=api.read_plan(parent)
    api.verify_sources(old["producer_sha256"]);api.verify_sources(old["input_sha256"])
    paths=[str(args.controls.resolve()/"controls"/str(m)/"status.json") for m in (1,2)]
    witnesses=checked_controls(paths,args.workspace_GiB)
    if (len(old["tasks"])!=1 or old["tasks"][0]["charge"]!=0 or old["tasks"][0]["mult"]!=1
            or args.minimum_GPU_MiB < (args.workspace_GiB+1)*1024):
        raise ValueError("One declared task and positive device headroom required")
    out.mkdir(parents=True,exist_ok=False)
    shutil.copytree(parent/"live_inputs",out/"live_inputs");shutil.copytree(parent/"controls",out/"controls")
    for stage in ("energy","thermal"):
        for name in ("tasks","done","claims","failures","workers","pilot","launches"):
            (out/stage/name).mkdir(parents=True)
    inputs={str(parent/"plan.json"):api.digest(parent/"plan.json")}
    failure=parent/"energy/pilot/status.json"
    if json.loads(failure.read_text())["status"]!="failed_preserved":
        raise ValueError("An observed parent resource failure is required")
    inputs[str(failure)]=api.digest(failure)
    reference=failure.parent/"energy_and_gradient.json"
    if reference.exists():inputs[str(reference)]=api.digest(reference)
    for task in old["tasks"]:
        path=out/"live_inputs"/task["key"]/"task.json";ready=json.loads(path.read_text())
        ready["input_xyz"]=str(path.parent/"input.xyz");ready["input_sha256"][ready["input_xyz"]]=api.digest(ready["input_xyz"])
        atomic_write_json(path,ready);inputs[str(path)]=api.digest(path);inputs.update(ready["input_sha256"])
    sources=dict(old["producer_sha256"],**witnesses)
    for path in (Path(__file__).resolve(),ADAPTER,ROOT/"scripts/pipeline/smd_bounded_execution.py",
                 ROOT/"scripts/pipeline/smd_memory_execution.py",ROOT/"scripts/pipeline/run_bounded_smd_recovery.py"):
        sources[str(path)]=api.digest(path)
    package=Path(importlib.util.find_spec("gpu4pyscf").origin).parent
    for path in (package/"df/df.py",package/"scf/hf.py",package/"solvent/smd.py"):
        sources[str(path)]=api.digest(path)
    plan=dict(old,raw_root=str(ROOT/"calcs/validation"/out.parent.name/out.name),input_sha256=inputs,producer_sha256=sources,
        source_failed_resource_pilot=str(parent/"energy/pilot/status.json"),execution_workspace_GiB=args.workspace_GiB,
        SMD_memory_qualification_controls=paths,minimum_large_GPU_MiB=args.minimum_GPU_MiB,
        maximum_memory_pilot_MiB=int(args.workspace_GiB*1024),
        purpose="Resource/numeric-continuity test at the exact starting geometry; no optimized state, G or synthesis admission",
        qualified_full_geometry_workflow=False,qualified_native_Hessian=False,prepared_at=time.time())
    atomic_write_json(out/"plan.json",plan)
    print(json.dumps(dict(status="isolated_footprint_pilot_prepared",tasks=1,workspace_GiB=args.workspace_GiB)))


def tensor_shapes(task):
    import cupy
    from ase.io import read
    from pyscf import gto
    from gpu4pyscf.df import df as gpu_df
    from redox.qm import dft
    atoms=read(task["input_xyz"])
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions)),basis="def2-tzvpd",charge=0,spin=0,unit="Angstrom",verbose=0)
    mf=dft._build_smd_mf(mol,"wb97m-v","acetonitrile","gpu",nlc="vv10")
    mf.with_solvent.build()
    original=gpu_df.empty_mapped
    class ShapeCaptured(Exception):pass
    found=[]
    def intercept(shape,*args,**kwargs):
        found.append(dict(shape=list(map(int,shape)),bytes=int(shape[0])*int(shape[1])*8))
        raise ShapeCaptured()
    gpu_df.empty_mapped=intercept
    try:
        omega=mf._numint.rsh_and_hybrid_coeff(mf.xc,spin=0)[0]
        for value in (None,float(omega)):
            with mf.with_df.range_coulomb(value) as dfobj:
                dfobj.use_gpu_memory=False
                try:dfobj.build()
                except ShapeCaptured:pass
                else:raise ValueError("DF shape capture did not intercept host tensor allocation")
    finally:
        gpu_df.empty_mapped=original
        del mf,mol
        gc.collect();cupy.get_default_memory_pool().free_all_blocks();cupy.cuda.runtime.deviceSynchronize()
    if len(found)!=2:raise ValueError("Both full/range-separated tensors must be inventoried")
    return found


def pilot(args):
    import psutil
    from redox.workflows.resources import verify_gate_owner,host_capacity,gpu_suitability
    verify_gate_owner()
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!="2" for k in
        ("OMP_NUM_THREADS","OMP_THREAD_LIMIT","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMEXPR_NUM_THREADS","XTB_THREADS")):
        raise ValueError("Footprint worker requires two disjoint cores and native caps=2")
    campaign=args.campaign.resolve();plan=api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"]);api.verify_sources(plan["input_sha256"])
    controls=checked_controls(plan["SMD_memory_qualification_controls"],plan["execution_workspace_GiB"])
    device=gpu_suitability(int(os.environ["CUDA_VISIBLE_DEVICES"]))
    if device["memory_MiB"]<plan["minimum_large_GPU_MiB"]:raise ValueError("Wrong resource-profile device class")
    driver=api.load_module("isolated_workspace_pilot",ROOT/"scripts/pipeline/refine_polymer_conformers.py")
    backend=driver.backend(campaign,"energy");task=backend.task_input(campaign,plan["tasks"][0])
    path=campaign/"footprint_qualification.json"
    if path.exists():raise FileExistsError("Preserve prior footprint attempt")
    report=dict(status="inventorying_resources",host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),
        GPU=os.environ["CUDA_VISIBLE_DEVICES"],started=time.time(),free_energy_computed=False,optimized_geometry_qualified=False,
        native_Hessian_qualified=False,synthesis_admission=False)
    atomic_write_json(path,report)
    adapter=api.load_module("tested_workspace_execution",ADAPTER)
    raw=Path(plan["raw_root"]);raw.mkdir(parents=True,exist_ok=False)
    try:
        with adapter.compact_smd_execution(raw/"shape_execution",budget_GiB=plan["execution_workspace_GiB"]):
            shapes=tensor_shapes(task)
        required=int(1.2*sum(r["bytes"] for r in shapes)+24*2**30)
        available=psutil.virtual_memory().available
        report.update(DF_tensor_shapes=shapes,host_required_bytes=required,available_host_RAM_bytes=available)
        if required>.75*available:
            raise RuntimeError("Screened DF tensor inventory exceeds conservative available host RAM; no SCF started")
        host_capacity(2,minimum_ram_gib=required/2**30)
        report["status"]="energy_gradient_resource_pilot_running";atomic_write_json(path,report)
        with adapter.compact_smd_execution(raw/"pilot_execution",budget_GiB=plan["execution_workspace_GiB"]):
            driver.calculate_energy(campaign,task,pilot=True)
        continuity=api.load_module("workspace_continuity_check",ROOT/"scripts/pipeline/run_bounded_smd_recovery.py")
        compared=continuity.compare_previous_pilot(plan,campaign/"energy/pilot")
        record=json.loads((campaign/"energy/pilot/status.json").read_text())
        report.update(status=("fixed_geometry_workspace_and_numeric_continuity_qualified" if compared["compared"]
            else "fixed_geometry_workspace_qualified_larger_reference_unavailable"),previous_reference=compared,
            observed_peak_GPU_MiB=record["observed_peak_GPU_MiB"],qualification_controls=controls)
    except BaseException as error:
        report.update(status="footprint_attempt_failed_preserved",error=repr(error));raise
    finally:
        report.update(finished=time.time(),host_peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        atomic_write_json(path,report)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("action",choices=("prepare","pilot"))
    parser.add_argument("--campaign",type=Path,required=True);parser.add_argument("--parent",type=Path)
    parser.add_argument("--controls",type=Path);parser.add_argument("--workspace-GiB",type=float,default=28)
    parser.add_argument("--minimum-GPU-MiB",type=int,default=31*1024);args=parser.parse_args()
    if args.action=="prepare":
        if not args.parent or not args.controls:parser.error("Preparation requires parent and controls")
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
        prepare(args)
    else:pilot(args)


if __name__=="__main__":main()
