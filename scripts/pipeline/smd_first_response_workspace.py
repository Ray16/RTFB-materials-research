"""Compose checked SMD execution with bounded PCM first-derivative rectangles.

This is unqualified until full original-versus-modified spin controls and an
exact large-system resource/derivative pilot pass. Installed files are unchanged.
"""
from contextlib import contextmanager
import hashlib
import inspect
from pathlib import Path
import time


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    import cupy
    from gpu4pyscf.solvent.hessian import pcm
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    source = ROOT / "scripts/pipeline/xc_workspace_execution.py"
    tensor_source = ROOT / "scripts/pipeline/pcm_first_derivative_blocks.py"
    underlay = api.load_module("first_response_XC_underlay", source)
    tensors = api.load_module("bounded_PCM_first_derivatives", tensor_source)
    directory = Path(directory)
    with underlay.compact_smd_execution(directory, budget_GiB=budget_GiB) as proof:
        path = directory / "first_derivatives.json"
        original = pcm.get_dD_dS
        if "dS = cupy.empty([3,n,n])" not in inspect.getsource(original):
            raise ValueError("Installed first-derivative primitive changed")
        statistics = dict(status="first_derivative_execution_active_pending_controls",
            started=time.time(), factories=0, primitive_calls=0, rectangles=0,
            maximum_output_bytes=0, maximum_primitive_points=0,
            maximum_primitive_tensor_bytes=0, tile_points=512,
            maximum_output_bound_bytes=512*1024**2,
            unchanged_pairwise_primitive=True, surface_grid_changed=False,
            precision_changed=False, installed_library_modified=False)
        digest = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        statistics["source_sha256"] = {str(p): digest(p) for p in
            (Path(__file__).resolve(), source, tensor_source, Path(inspect.getfile(original)))}

        def bounded(surface, with_S=True, with_D=False, stream=None):
            statistics["factories"] += 1
            primitive = lambda local, **kw: original(local, stream=stream, **kw)
            blocks = tensors.FirstDerivativeBlocks(surface, primitive, cupy, statistics)
            return (tensors.FirstDerivativeTensor(blocks, 0) if with_D else None,
                    tensors.FirstDerivativeTensor(blocks, 1))

        pcm.get_dD_dS = bounded
        proof["source_sha256"].update(statistics["source_sha256"])
        proof.update(PCM_first_derivative_blocks=True, PCM_first_derivative_equations_changed=False)
        atomic_write_json(path, statistics)
        try:
            yield proof
            statistics["status"] = "first_derivative_execution_finished_requires_controls_and_physical_review"
        except BaseException as error:
            statistics.update(status="first_derivative_execution_failed_preserved", error=repr(error))
            raise
        finally:
            pcm.get_dD_dS = original
            statistics["finished"] = time.time()
            atomic_write_json(path, statistics)
