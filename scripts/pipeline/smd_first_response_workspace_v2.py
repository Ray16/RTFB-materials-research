"""Bound first derivatives and apply Hessian R operators without dense R.

The dense reference equations are retained through the installed solvent
operator. Independent full spin controls and a large response test are required.
"""
from contextlib import contextmanager
import hashlib
import inspect
from pathlib import Path
import textwrap
import time


def transform(name, body):
    if name == "get_dqsym_dx_fix_vgrids":
        old = """        R = -f_epsilon * (cupy.eye(DA.shape[0]) - 1.0/(2.0*PI)*DA)
        dqdx_fix_Vq += -contract('ij,Adj->Adi', R.T, invKT_dKdxT_dot_invKT_V)"""
        if body.count(old) != 1:
            raise ValueError("Installed fixed-potential response template changed")
        return body.replace(old, """        del DA, AS
        dqdx_fix_Vq -= _apply_R_batch(pcmobj, invKT_dKdxT_dot_invKT_V, True)""", 1)
    if name == "get_dqsym_dx_fix_K_R":
        old = """        A = pcmobj._intermediates['A']
        D = pcmobj._intermediates['D']
        DA = D * A
        R = -f_epsilon * (cupy.eye(DA.shape[0]) - 1.0/(2.0*PI)*DA)
        R_dVdx = contract('ij,Adj->Adi', R, dV_on_charge_dx)"""
        other = "        RT_K_1T_dVdx = contract('ij,Adj->Adi', R.T, K_1T_dVdx)"
        if body.count(old) != 1 or body.count(other) != 1:
            raise ValueError("Installed fixed-operator response template changed")
        return body.replace(old, "        R_dVdx = _apply_R_batch(pcmobj, dV_on_charge_dx, False)").replace(
            other, "        RT_K_1T_dVdx = _apply_R_batch(pcmobj, K_1T_dVdx, True)")
    if name == "analytical_hess_solver":
        old = """        DA = D*A
        R = -f_epsilon * (cupy.eye(DA.shape[0]) - 1.0/(2.0*PI)*DA)
        dvK_1R = -einsum_Adi_ij_Adj_inverseK(vK_1_dot_dKdx, pcmobj) @ R + VK_1_dot_dRdx"""
        if body.count(old) != 2:
            raise ValueError("Installed analytic-Hessian R template changed")
        return body.replace(old, """        dvK_1R = -_apply_R_batch(pcmobj,
            einsum_Adi_ij_Adj_inverseK(vK_1_dot_dKdx, pcmobj), True) + VK_1_dot_dRdx""", 1)
    raise ValueError("Unsupported PCM transformation")


def apply_R_batch(solvent, values, transpose):
    shape = values.shape
    if len(shape) != 3 or shape[1] != 3:
        raise ValueError("Expected atom/Cartesian/surface response tensor")
    return solvent.left_multiply_R(values.reshape(-1, shape[-1]).T,
        R_transpose=transpose).T.reshape(shape)


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    from gpu4pyscf.solvent.hessian import pcm
    from gpu4pyscf.solvent import smd
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    source = ROOT / "scripts/pipeline/smd_first_response_workspace.py"
    base = api.load_module("PCM_first_derivative_execution_underlay", source)
    with base.compact_smd_execution(directory, budget_GiB=budget_GiB) as proof:
        names = ("get_dqsym_dx_fix_vgrids", "get_dqsym_dx_fix_K_R", "analytical_hess_solver")
        originals = {name: getattr(pcm, name) for name in names}
        stats = dict(status="contracted_PCM_Hessian_R_active_pending_controls", started=time.time(),
            transformed_calls={name: 0 for name in names}, R_applications=0,
            equations_changed=False, precision_changed=False, surface_grid_changed=False)

        def apply(solvent, value, transpose):
            stats["R_applications"] += 1
            return apply_R_batch(solvent, value, transpose)

        namespace = dict(vars(pcm), _apply_R_batch=apply)
        replacements = {}
        for name, original in originals.items():
            before = inspect.getsource(original)
            after = transform(name, before)
            path = Path(directory) / ("contracted_"+name+".py")
            header = "\n".join(Path(pcm.__file__).read_text().splitlines()[:13])+"\n\n"
            path.write_text(header+after)
            local = {}
            exec(compile(textwrap.dedent(after), str(path), "exec"), namespace, local)
            replacement = local[name]

            def wrapped(obj, *args, _name=name, _old=original, _new=replacement, **kwargs):
                if not isinstance(obj, smd.SMD):
                    return _old(obj, *args, **kwargs)
                stats["transformed_calls"][_name] += 1
                return _new(obj, *args, **kwargs)

            replacements[name] = wrapped
            setattr(pcm, name, wrapped)
        namespace.update(replacements)
        own = Path(__file__).resolve()
        sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        stats["source_sha256"] = {str(p): sha(p) for p in (own, source, Path(pcm.__file__))}
        proof["source_sha256"].update(stats["source_sha256"])
        proof.update(PCM_Hessian_R_contracted=True, PCM_Hessian_R_equations_changed=False)
        try:
            yield proof
            stats["status"] = "contracted_PCM_Hessian_R_finished_requires_controls_and_physical_review"
        except BaseException as error:
            stats.update(status="contracted_PCM_Hessian_R_failed_preserved", error=repr(error))
            raise
        finally:
            for name, original in originals.items():
                setattr(pcm, name, original)
            stats["finished"] = time.time()
            atomic_write_json(Path(directory) / "contracted_R.json", stats)
