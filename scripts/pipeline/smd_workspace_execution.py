"""Budget-aware SMD execution with persistent-solvent-first allocation.

Reuse the tested contraction/buffer adapter. Only allocation order, workspace
estimates and durable SCF logging change; no physical or numerical grid changes.
"""
from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path
import time


def molecular_identity(mf):
    solvent = getattr(mf, "with_solvent", None)
    return dict(atoms=list(mf.mol.elements), coordinates_Bohr=mf.mol.atom_coords().tolist(),
        basis=repr(mf.mol.basis), charge=int(mf.mol.charge), mult=int(mf.mol.spin)+1,
        xc=str(getattr(mf,"xc",None)), disp=str(getattr(mf,"disp",None)), nlc=str(getattr(mf,"nlc",None)),
        solvent=getattr(solvent,"solvent",None), init_guess=str(getattr(mf,"init_guess",None)))


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    from redox.qm import dft
    from redox.core.io import atomic_write_json
    directory = Path(directory)
    directory.mkdir(parents=True,exist_ok=False)
    source = Path(__file__).with_name("smd_bounded_execution.py")
    spec = importlib.util.spec_from_file_location("bounded_workspace_underlay",source)
    bounded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bounded)
    sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    report = dict(status="workspace_execution_active_pending_qualification",started=time.time(),budget_GiB=budget_GiB,
        prebuilt_solvent_calls=0,SCF_calls=[],source_sha256={str(p):sha(p) for p in (Path(__file__).resolve(),source,Path(dft.__file__))},
        scientific_model_changed=False,surface_grid_changed=False,precision_changed=False,
        cached_SCF_data_are_qualified_results=False)
    try:
        with bounded.compact_smd_execution(directory/"bounded",budget_GiB=budget_GiB):
            original = dft._kernel_robust
            def kernel(mf,dm0=None):
                index = len(report["SCF_calls"])
                folder = directory/"SCF_checkpoints"/f"{index:03d}"
                folder.mkdir(parents=True,exist_ok=False)
                identity = molecular_identity(mf)
                progress = dict(identity=identity,qualified_result=False,SCF_converged=False,phase="starting",started=time.time())
                atomic_write_json(folder/"progress.json",progress)
                mf.chkfile = str(folder/"SCF.h5")
                previous_callback = getattr(mf,"callback",None)
                def callback(environment):
                    if callable(previous_callback):
                        previous_callback(environment)
                    progress.update(phase="SCF_iteration",cycle=int(environment["cycle"])+1,
                        E_Ha=float(environment["e_tot"]),SCF_converged=False,updated=time.time())
                    atomic_write_json(folder/"progress.json",progress)
                mf.callback = callback
                report["SCF_calls"].append(dict(path=str(folder),identity=identity))
                solvent = getattr(mf,"with_solvent",None)
                if solvent is not None and not getattr(solvent,"_intermediates",None):
                    solvent.build()
                    report["prebuilt_solvent_calls"] += 1
                atomic_write_json(directory/"status.json",report)
                try:
                    result = original(mf,dm0=dm0)
                    progress.update(phase="SCF_returned",SCF_converged=bool(result[1]),E_Ha=float(result[0]),
                        conv_tol=float(mf.conv_tol),finished=time.time())
                    if Path(mf.chkfile).exists():
                        progress["checkpoint_sha256"] = sha(mf.chkfile)
                    return result
                except BaseException as error:
                    progress.update(phase="SCF_failed_checkpoint_unqualified",error=repr(error),finished=time.time())
                    raise
                finally:
                    mf.callback = previous_callback
                    atomic_write_json(folder/"progress.json",progress)
            dft._kernel_robust = kernel
            try:
                yield report
            finally:
                dft._kernel_robust = original
        report["status"] = "workspace_execution_finished_requires_control_and_physical_review"
    except BaseException as error:
        report.update(status="workspace_execution_failed_preserved",error=repr(error))
        raise
    finally:
        report["finished"] = time.time()
        atomic_write_json(directory/"status.json",report)
