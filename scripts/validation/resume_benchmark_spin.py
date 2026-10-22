#!/usr/bin/env python
"""Resume the saved spin-helper geometry, preserving normal energy and spin gates."""
import argparse
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

PREVIOUS = ROOT / "results/validation/publication_benchmark/spin_completion_2026-10-21"
OUT = ROOT / "results/validation/publication_benchmark/spin_resume_2026-10-22"
GID = "bm_n_n_bis_2_2_methoxyethoxy_ethyl_pyromellitic_diimide"
STATE = "red2"
SCRIPT = Path(__file__).resolve()
CAPS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")


def geometry_evidence(folder):
    """A saved callback alone never establishes optimizer convergence."""
    import numpy as np
    from ase.io import read
    returned = json.loads((folder / "optimizer_returned.json").read_text())
    last = json.loads((folder / "latest.json").read_text())
    if returned.get("optimizer_converged") is not True or last.get("SCF_converged") is not True:
        raise ValueError("A real converged optimizer return and SCF are required")
    a, b = read(folder / "optimizer_returned.xyz"), read(folder / "last_evaluated.xyz")
    if a.get_chemical_symbols() != b.get_chemical_symbols() or not np.allclose(a.positions, b.positions, rtol=0, atol=2e-7):
        raise ValueError("Optimizer return does not match the converged final SCF evaluation")
    return dict(optimizer_converged=True, SCF_converged=True,
                source_sha256={str(folder / name): api.digest(folder / name) for name in
                               ("optimizer_returned.json", "optimizer_returned.xyz", "latest.json", "last_evaluated.xyz")})


def prepare():
    from redox.validation.spins import snapshot
    from redox.validation.spin_evidence import spin_comparison_running
    if spin_comparison_running(GID):
        raise RuntimeError("An existing spin writer is active")
    previous = json.loads((PREVIOUS / "status.json").read_text())
    if previous.get("status") != "failed_preserved" or "MemoryAllocation" not in previous["result"].get("error", ""):
        raise ValueError("Expected preserved GPU-memory failure")
    evidence = geometry_evidence(PREVIOUS / "geometry_checkpoints/optimizer_0")
    source = snapshot(GID, STATE)
    if source != json.loads((PREVIOUS / "plan.json").read_text())["chosen_identity"]:
        raise ValueError("Chosen physical state changed")
    paths = [Path(p) for p in json.loads((PREVIOUS / "plan.json").read_text())["producer_sha256"]]
    paths += [SCRIPT, PREVIOUS / "status.json", PREVIOUS / "plan.json"]
    inputs = evidence["source_sha256"]
    OUT.mkdir()
    atomic_write_json(OUT / "plan.json", dict(id=GID, state=STATE, chosen_identity=source,
        producer_sha256={str(p): api.digest(p) for p in paths}, input_sha256=inputs,
        optimizer_evidence=evidence, minimum_GPU_memory_MiB=81920,
        numerical_settings="Original active SP, open-shell guess sweep and source-addressed spin/native gates",
        role="Exact-parent method control; no synthesis admission", maximum_actual_executions=1))


@contextmanager
def reuse_geometry(checkpoint, helper_geometry):
    """Score a genuine optimizer result without retaining its large SCF objects."""
    from redox.qm import dft
    previous = dft.dft_smd
    consumed = False
    evidence = geometry_evidence(checkpoint)

    def resume(xyz, charge, mult, **kwargs):
        nonlocal consumed
        if kwargs.get("opt_out") == helper_geometry and not consumed:
            if charge != -2 or mult != 3 or kwargs != dict(backend="gpu", opt_out=helper_geometry, do_freq=False):
                raise ValueError("Saved optimizer geometry may only resume the declared helper call")
            if Path(xyz) != ROOT / "calcs/dft" / GID / STATE / "opt.xyz":
                raise ValueError("Unexpected chosen input for saved spin helper")
            consumed = True
            helper_geometry.write_bytes((checkpoint / "optimizer_returned.xyz").read_bytes())
            result = previous(helper_geometry, charge, mult, backend="gpu", do_opt=False, do_freq=False)
            result.update(optimized=True, geometry_optimization_converged=True,
                          optimization_SCF_converged=True, opt_gradient_grid_response=True,
                          opt_out=str(helper_geometry), optimization_checkpoint=evidence,
                          optimization_reuse="Genuine converged optimizer and matching SCF; ordinary native QC still required")
            return result
        return previous(xyz, charge, mult, **kwargs)

    dft.dft_smd = resume
    try:
        yield
    finally:
        dft.dft_smd = previous


def worker():
    from redox.workflows.resources import verify_gate_owner, gpu_suitability
    from redox.validation.spins import snapshot, check
    from redox.validation.spin_evidence import state_spin_evidence
    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != "2" for k in CAPS):
        raise RuntimeError("Two disjoint CPUs and native caps=2 required")
    if gpu_suitability(int(os.environ["CUDA_VISIBLE_DEVICES"]))["memory_MiB"] < 81920:
        raise RuntimeError("The failed 32-GB run requires a larger device")
    plan = json.loads((OUT / "plan.json").read_text())
    api.verify_sources(plan["producer_sha256"]); api.verify_sources(plan["input_sha256"])
    if snapshot(GID, STATE) != plan["chosen_identity"]:
        raise ValueError("Chosen source changed")
    original = api.load_module("checkpointed_spin_original", ROOT / "scripts/validation/complete_benchmark_spin.py")
    original.OUT = OUT
    checkpoint = PREVIOUS / "geometry_checkpoints/optimizer_0"
    helper = ROOT / "calcs/spincheck" / GID / "red2_m3/opt.xyz"
    status = dict(status="computing", host=socket.gethostname(), pid=os.getpid(), started=time.time(),
                  process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                  gpu=os.environ["CUDA_VISIBLE_DEVICES"], cpus=sorted(os.sched_getaffinity(0)))
    atomic_write_json(OUT / "status.json", status)
    try:
        with original.checkpoint_optimizers(), reuse_geometry(checkpoint, helper):
            result = check(GID, STATE)
        evidence = state_spin_evidence(GID, STATE)
        status.update(result=result, publication_spin_evidence=evidence,
                      status="spin_admission_verified" if evidence["spin_primary_admissible"] else "requires_physical_review")
        if result["status"] in ("failed", "source_stale"):
            raise RuntimeError("Spin continuation failed; preserved for review")
    except BaseException as error:
        status.update(status="failed_preserved", error=repr(error))
        raise
    finally:
        status["finished"] = time.time()
        atomic_write_json(OUT / "status.json", status)


def coordinate(max_hours):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    with (OUT / "coordinator.lock").open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = dict(host=socket.gethostname(), pid=os.getpid(), started=time.time(),
                     process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                     deadline=time.time()+max_hours*3600, CPU_affinity=sorted(os.sched_getaffinity(0)))
        while time.time() < owner['deadline']:
            plan=json.loads((OUT/'plan.json').read_text())
            api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
            if (OUT/'status.json').exists():
                atomic_write_json(OUT/'coordinator.json',dict(owner,status='worker_started_or_terminal_do_not_duplicate'))
                return
            available=free_slots()
            for host in sorted({h for h,_ in available}):
                p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,
                    'nvidia-smi --query-gpu=index,memory.total --format=csv,noheader,nounits'],capture_output=True,text=True,timeout=20)
                memory={int(r.split(',')[0]):float(r.split(',')[1]) for r in p.stdout.splitlines()} if p.returncode==0 else {}
                eligible=[g for h,g in available if h==host and memory.get(g,0)>=81920]
                if not eligible:continue
                preflight=[sys.executable,str(ROOT.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'),host,
                    '--python',sys.executable,'--module','pyscf','--module','gpu4pyscf','--module','ase','--module','rdkit',
                    '--path',str(OUT/'plan.json'),'--min-available-ram-mib','32768','--timeout','40']
                p=subprocess.run(preflight,capture_output=True,text=True,timeout=50)
                (OUT/f'{host}_preflight.json').write_text(p.stdout)
                if p.returncode:continue
                destination=OUT/'execution'/str(time.time_ns())
                argv=['env','PYTHONPATH='+str(ROOT/'src'),*[k+'=1' for k in CAPS],sys.executable,
                    '-m','redox.workflows.gpu','--gpu',str(eligible[0]),'--threads','2','--out',str(destination),
                    '--',sys.executable,str(SCRIPT),'worker']
                log=OUT/f'launch_{destination.name}.log'
                command='nohup '+shlex.join(argv)+' >'+shlex.quote(str(log))+' 2>&1 </dev/null & echo $!'
                p=subprocess.run(['ssh','-o','BatchMode=yes',host,command],capture_output=True,text=True,check=True,timeout=20)
                launch=dict(host=host,pid=p.stdout.strip(),argv=argv,execution=str(destination),started=time.time())
                atomic_write_json(OUT/f'launch_{destination.name}.json',launch)
                while time.time()<owner['deadline']:
                    if (OUT/'status.json').exists():
                        atomic_write_json(OUT/'coordinator.json',dict(owner,status='worker_started',launch=launch))
                        return
                    execution=destination/'execution.json'
                    if execution.exists() and json.loads(execution.read_text()).get('status') in ('failed','process_completed'):
                        receipt=json.loads(execution.read_text())
                        if 'occupied or reserved; no calculation started' not in receipt.get('error',''):
                            atomic_write_json(OUT/'coordinator.json',dict(owner,status='execution_failed_review_required',launch=launch,execution=receipt))
                            return
                        break
                    time.sleep(15)
                break
            atomic_write_json(OUT/'coordinator.json',dict(owner,status='waiting_all_suitable_GPUs_occupied_or_contention',updated=time.time()))
            time.sleep(30)
        atomic_write_json(OUT/'coordinator.json',dict(owner,status='finite_resource_window_ended_outputs_preserved'))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','worker','coordinate'))
    parser.add_argument('--max-hours',type=float,default=24)
    args=parser.parse_args()
    if not 0<args.max_hours<=48:parser.error('Finite window up to 48 hours required')
    if args.action=='prepare':prepare()
    elif args.action=='worker':worker()
    else:coordinate(args.max_hours)
