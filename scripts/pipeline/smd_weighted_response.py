"""Apply diagonal-weighted PCM matrices without surface-square temporaries.

FP64 row tiles preserve the installed SMD equations and full surface. This
adapter requires original-versus-adapted full spin/Hessian controls and a
separate resource/derivative check on each newly qualified large model.
"""
from contextlib import contextmanager
import inspect
from pathlib import Path
import textwrap
import time


class ScaledMatrix:
    """Lazy diag(left) @ matrix @ diag(right), with bounded row products."""

    def __init__(self, matrix, xp, left=None, right=None, tile=256, metrics=None):
        if len(matrix.shape) != 2 or matrix.dtype != xp.float64 or tile < 1:
            raise ValueError("A two-dimensional FP64 matrix and positive tile are required")
        for weight, length in ((left, matrix.shape[0]), (right, matrix.shape[1])):
            if weight is not None and (weight.shape != (length,) or weight.dtype != xp.float64):
                raise ValueError("FP64 diagonal weights must match the matrix axes")
        self.matrix, self.xp, self.left, self.right = matrix, xp, left, right
        self.tile, self.metrics = tile, metrics
        self.shape, self.dtype, self.ndim = matrix.shape, matrix.dtype, 2

    @property
    def T(self):
        return ScaledMatrix(self.matrix.T, self.xp, self.right, self.left,
                            self.tile, self.metrics)

    def __matmul__(self, value):
        xp = self.xp
        if value.ndim not in (1, 2) or value.shape[0] != self.shape[1] or value.dtype != xp.float64:
            raise ValueError("FP64 vector or column batch with matching surface axis required")
        vector = value.ndim == 1
        batch = value[:, None] if vector else value
        if self.right is not None:
            batch = batch * self.right[:, None]
        result = xp.empty((self.shape[0], batch.shape[1]), dtype=xp.float64)
        for start in range(0, self.shape[0], self.tile):
            end = min(start + self.tile, self.shape[0])
            block = xp.ascontiguousarray(self.matrix[start:end, :])
            result[start:end] = block @ batch
            if self.left is not None:
                result[start:end] *= self.left[start:end, None]
            if self.metrics is not None:
                self.metrics['matrix_tiles'] += 1
                self.metrics['maximum_matrix_tile_bytes'] = max(
                    self.metrics['maximum_matrix_tile_bytes'], block.nbytes)
            del block
        return result[:, 0] if vector else result


def response_contract(subscripts, first, second, *, xp, fallback, metrics):
    """Bound the two matrix/Cartesian-batch contractions used by PCM."""
    if subscripts == 'ij,Adj->Adi':
        matrix, response = first, second
        transpose = False
    elif subscripts == 'Adi,ij->Adj':
        response, matrix = first, second
        transpose = True
    else:
        return fallback(subscripts, first, second)
    if len(response.shape) != 3 or response.shape[1] != 3:
        raise ValueError("Expected atom/Cartesian/surface response tensor")
    if not isinstance(matrix, ScaledMatrix):
        matrix = ScaledMatrix(matrix, xp, metrics=metrics)
    if transpose:
        matrix = matrix.T
    metrics['bounded_response_contracts'] += 1
    result = matrix @ response.reshape(-1, response.shape[-1]).T
    return result.T.reshape(response.shape[0], 3, matrix.shape[0])


def transformed_source(body, name='get_dqsym_dx_fix_vgrids'):
    # Both branches exist upstream; only the IEFPCM/SMD branch is qualified.
    changes = (('        DA = D*A', '        DA = _scaled_matrix(D, right=A)'),
               ("        AS = (A * S.T).T # It's just diag(A) @ S",
                '        AS = _scaled_matrix(S, left=A)'))
    for old, new in changes:
        expected = 3 if name == 'analytical_hess_solver' and 'DA = D*A' in old else 2
        if name not in ('get_dqsym_dx_fix_vgrids', 'analytical_hess_solver') or body.count(old) != expected:
            raise ValueError('Installed SMD weighted-response template changed')
        body = body.replace(old, new, 1)
    return body


@contextmanager
def weighted_response(directory):
    """Compose around the v2 first response or the streamed-curvature profile."""
    import cupy as cp
    from gpu4pyscf.solvent.hessian import pcm
    from gpu4pyscf.solvent import smd
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    names = ('get_dqsym_dx_fix_vgrids', 'get_dqsym_dx_fix_K_R', 'analytical_hess_solver')
    originals = {name: getattr(pcm, name) for name in names}
    clones = [getattr(function, '__kwdefaults__', {}).get('_new') for function in originals.values()]
    if any(clone is None for clone in clones):
        raise ValueError('The checked first-response underlay must be active')
    scopes = {id(clone.__globals__): clone.__globals__ for clone in clones}.values()
    metrics = dict(status='weighted_response_active_unqualified', started=time.time(),
        weighted_matrices=0, bounded_response_contracts=0, matrix_tiles=0,
        transformed_calls={name: 0 for name in (names[0], names[2])},
        maximum_matrix_tile_bytes=0, tile_rows=256, precision='float64',
        equations_changed=False, surface_grid_changed=False, atoms_omitted=False,
        source_sha256={str(Path(__file__).resolve()): api.digest(__file__),
                       str(Path(pcm.__file__)): api.digest(pcm.__file__)})
    fallback = clones[0].__globals__['contract']

    def contract(subscripts, first, second, **kwargs):
        if kwargs:
            return fallback(subscripts, first, second, **kwargs)
        return response_contract(subscripts, first, second, xp=cp,
                                 fallback=fallback, metrics=metrics)

    def scaled(matrix, left=None, right=None):
        metrics['weighted_matrices'] += 1
        return ScaledMatrix(matrix, cp, left, right, metrics=metrics)

    header = '\n'.join(Path(pcm.__file__).read_text().splitlines()[:13]) + '\n\n'
    namespace = dict(clones[0].__globals__, contract=contract, _scaled_matrix=scaled)
    replacements = {}
    for name, clone in ((names[0], clones[0]), (names[2], clones[2])):
        source = transformed_source(inspect.getsource(clone), name)
        path = directory / ('weighted_' + name + '.py')
        path.write_text(header + source)
        local = {}
        exec(compile(textwrap.dedent(source), str(path), 'exec'), namespace, local)
        replacement = local[name]

        def wrapped(obj, *args, _new=replacement, _name=name, **kwargs):
            if not isinstance(obj, smd.SMD) or obj.method.upper() != 'SMD':
                return originals[_name](obj, *args, **kwargs)
            metrics['transformed_calls'][_name] += 1
            return _new(obj, *args, **kwargs)

        replacements[name] = wrapped
    previous = [(scope, {key: scope[key] for key in ('contract', *replacements)}) for scope in scopes]
    for scope, _ in previous:
        scope.update(contract=contract, **replacements)
    namespace.update(replacements)
    for name, function in replacements.items():
        setattr(pcm, name, function)
    try:
        yield metrics
        metrics['status'] = 'weighted_response_finished_requires_controls_and_physical_review'
    except BaseException as error:
        metrics.update(status='weighted_response_failed_preserved', error=repr(error))
        raise
    finally:
        for name in replacements:
            setattr(pcm, name, originals[name])
        for scope, before in previous:
            scope.update(before)
        metrics['finished'] = time.time()
        atomic_write_json(directory / 'weighted_response.json', metrics)


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    from redox.core.common import ROOT
    from redox.workflows import polymer_campaign as api
    underlay = api.load_module('weighted_response_v2_underlay',
        ROOT / 'scripts/pipeline/smd_first_response_workspace_v2.py')
    with underlay.compact_smd_execution(directory, budget_GiB=budget_GiB) as proof:
        with weighted_response(directory) as metrics:
            proof['source_sha256'].update(metrics['source_sha256'])
            proof.update(weighted_PCM_response=True, dense_weighted_surface_matrices=False)
            yield proof
