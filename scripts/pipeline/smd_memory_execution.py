"""Run-local temporary-buffer reduction for the qualified SMD equations.

The installed library is never edited. Qualification must compare full energies,
forces and Hessians before this execution adapter is used for scientific results.
No surface point, descriptor, functional, basis or numerical precision is changed.
"""
from contextlib import contextmanager
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time


def compact_build_source(source):
    old = """        DA = D*A
        DAS = cupy.dot(DA, S)
        K = S - f_epsilon/(2.0*np.pi) * DAS
        K_LU, K_LU_pivot = lu_factor(K, overwrite_a = True, check_finite = False)
        K = None"""
    new = """        DA = D*A
        DAS = cupy.dot(DA, S)
        del DA
        DAS *= -f_epsilon/(2.0*np.pi)
        DAS += S
        K_LU, K_LU_pivot = lu_factor(DAS, overwrite_a = True, check_finite = False)
        DAS = None"""
    if source.count(old) != 1:
        raise ValueError("Installed SMD build differs from the reviewed buffer-lifetime pattern")
    return source.replace(old, new)


def release_geometry_source(source):
    old = "        atoms = opt_atoms\n"
    if source.count(old) != 1:
        raise ValueError("Composite driver differs from the reviewed geometry/SP transition")
    return source.replace(old, old + "        del gradient, scanner, mf_opt\n        import gc\n        gc.collect()\n        _free_gpu_pool()\n")


@contextmanager
def compact_smd_execution(directory):
    from redox.workflows.resources import verify_gate_owner
    from redox.core.io import atomic_write_json
    from redox.qm.native_thermal import installed_versions
    from redox.qm import dft
    from gpu4pyscf.solvent import smd

    verify_gate_owner()
    if installed_versions() != {"pyscf": "2.14.0", "gpu4pyscf": "1.8.1"}:
        raise ValueError("Memory adapter requires the reviewed library versions")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    original_build, original_dft = smd.SMD.build, dft.dft_smd
    build_source = inspect.getsource(original_build)
    dft_source = inspect.getsource(original_dft)
    compact = compact_build_source(build_source)
    released = release_geometry_source(dft_source)
    namespace = {}
    exec(compile(textwrap.dedent(compact), str(directory / "compact_build.py"), "exec"), original_build.__globals__, namespace)
    optimized_build = namespace["build"]
    namespace = {}
    exec(compile(textwrap.dedent(released), str(directory / "released_geometry_driver.py"), "exec"), original_dft.__globals__, namespace)
    optimized_dft = namespace["dft_smd"]
    for name, content in (("original_build.py", build_source), ("compact_build.py", compact),
                          ("original_composite.py", dft_source), ("released_geometry_driver.py", released)):
        license_header = "\n".join(Path(smd.__file__).read_text().splitlines()[:13]) + "\n\n" if name.endswith("build.py") else ""
        (directory / name).write_text(license_header + content)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    proof = dict(status="execution_adapter_active_unqualified_until_controls_pass", started=time.time(), builds=[],
        source_sha256={str(p): sha(p) for p in (Path(__file__).resolve(), Path(smd.__file__), Path(dft.__file__))},
        scientific_model_changed=False, surface_grid_changed=False, precision_changed=False,
        transformations=["Release DA scratch before K assembly", "Reuse DAS storage for K and release its alias after LU", "Release unused geometry scanner/SCF before the larger SP"])
    atomic_write_json(directory / "status.json", proof)

    def build(solvent, ng=None):
        started = time.time()
        result = optimized_build(solvent, ng)
        proof["builds"].append(dict(atoms=solvent.mol.natm, surface_points=len(solvent.surface["charge_exp"]),
            persistent_matrix_bytes=sum(value.nbytes for value in solvent._intermediates.values() if hasattr(value, "nbytes")),
            elapsed_seconds=time.time() - started))
        atomic_write_json(directory / "status.json", proof)
        return result

    smd.SMD.build, dft.dft_smd = build, optimized_dft
    try:
        yield proof
    except BaseException as error:
        proof.update(status="failed_preserving_execution_evidence", error=repr(error))
        raise
    else:
        proof["status"] = "execution_adapter_finished_requires_control_and_physical_review"
    finally:
        smd.SMD.build, dft.dft_smd = original_build, original_dft
        proof["finished"] = time.time()
        atomic_write_json(directory / "status.json", proof)
