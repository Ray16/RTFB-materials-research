"""Run-local solvent contraction and workspace limits; unchanged model equations.

Small spin controls and a separate larger-system resource pilot must pass before
this execution variant can support an admitted scientific result.
"""
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time


def budget_available(total, free, pooled_free, budget):
    return max(0, int(budget - (total-free) + pooled_free))


def contracted_R(xp, D, A, vector, factor, transpose=False):
    """Apply D diag(A) to a vector/batch without constructing the dense product."""
    weights = A.reshape((len(A),) + (1,)*(vector.ndim-1))
    product = weights * xp.dot(D.T, vector) if transpose else xp.dot(D, weights*vector)
    return -factor*(vector-product/(2*math.pi))


@contextmanager
def compact_smd_execution(directory, budget_GiB=80):
    import cupy
    from gpu4pyscf.solvent import smd
    from gpu4pyscf.lib import cupy_helper
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    if not math.isfinite(budget_GiB) or budget_GiB <= 0:
        raise ValueError("Positive finite execution budget required")
    source = Path(__file__).with_name("smd_memory_execution.py")
    spec = importlib.util.spec_from_file_location("previous_qualified_SMD_build_adapter", source)
    previous = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(previous)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    own = Path(__file__).resolve()
    sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    proof = dict(status="bounded_execution_active_pending_qualification", started=time.time(), budget_GiB=budget_GiB,
        memory_queries=0, contracted_R_calls=0, source_sha256={str(p):sha(p) for p in (own,source)},
        installed_library_modified=False, surface_grid_changed=False, precision_changed=False,
        scientific_equivalence_requires_controls=True)
    with previous.compact_smd_execution(directory / "build_execution"):
        original_available = cupy_helper.get_avail_mem
        original_R = smd.SMD.left_multiply_R
        had_own_R = "left_multiply_R" in smd.SMD.__dict__
        def available(exclude_memory_pool=False):
            free, total = cupy.cuda.runtime.memGetInfo()
            pooled = 0 if exclude_memory_pool else cupy.get_default_memory_pool().free_bytes()
            proof["memory_queries"] += 1
            return min(original_available(exclude_memory_pool), budget_available(total, free, pooled, budget_GiB*2**30))
        def multiply(solvent, vector, R_transpose=False):
            if solvent.if_method_in_CPCM_category or getattr(solvent,"lowmem_intermediate_storage",False):
                return original_R(solvent,vector,R_transpose)
            proof["contracted_R_calls"] += 1
            data = solvent._intermediates
            return contracted_R(cupy,data["D"],data["A"],vector,data["f_epsilon"],R_transpose)
        bindings = []
        for name,module in list(sys.modules.items()):
            if name.startswith("gpu4pyscf") and module is not None:
                for key,value in list(vars(module).items()):
                    if value is original_available:
                        setattr(module,key,available)
                        bindings.append(name+":"+key)
        proof["bounded_memory_bindings"] = bindings
        smd.SMD.left_multiply_R = multiply
        atomic_write_json(directory / "status.json",proof)
        try:
            yield proof
        except BaseException as error:
            proof.update(status="bounded_execution_failed_preserved",error=repr(error))
            raise
        else:
            proof["status"] = "bounded_execution_finished_requires_control_and_physical_review"
        finally:
            if had_own_R:
                smd.SMD.left_multiply_R = original_R
            else:
                del smd.SMD.left_multiply_R
            for name,module in list(sys.modules.items()):
                if name.startswith("gpu4pyscf") and module is not None:
                    for key,value in list(vars(module).items()):
                        if value is available:
                            setattr(module,key,original_available)
            proof["finished"] = time.time()
            atomic_write_json(directory / "status.json",proof)
