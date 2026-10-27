"""Explicit relocation of immutable native-certificate provenance paths.

The original reader still checks every physical and numerical field. This layer
allows only a declared old checkout prefix to change; certificate contents and
all other metadata must agree exactly. It never opens the old checkout.
"""
from contextlib import contextmanager
import hashlib
import inspect
from pathlib import Path, PurePosixPath
import threading


def controls_equal(stored, expected, root, prior_root, receipts):
    if not isinstance(stored, list) or len(stored) != len(expected):
        return False
    root = Path(root).resolve()
    prior = PurePosixPath(str(prior_root))
    if not prior.is_absolute() or '..' in prior.parts:
        raise ValueError('Declare an absolute original checkout prefix without parent traversal')
    for old, current in zip(stored, expected):
        if not isinstance(old, dict) or not isinstance(current, dict):
            return False
        try:
            target = Path(current['path']).resolve()
            relative = target.relative_to(root)
        except (KeyError, ValueError):
            return False
        if (target == root or old.get('path') != str(prior / relative.as_posix())
                or {k:v for k,v in old.items() if k != 'path'} != {k:v for k,v in current.items() if k != 'path'}):
            return False
        try:
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
        except OSError:
            return False
        if actual != current.get('sha256'):
            return False
        receipts[str(relative)] = dict(original_path=old['path'], relocated_path=str(target), sha256=actual,
                                       metadata_unchanged_except_declared_prefix=True)
    return True


@contextmanager
def relocated_native_reader(prior_root):
    if threading.current_thread() is not threading.main_thread() or threading.active_count() != 1:
        raise RuntimeError('Install relocation binding before starting Python workers')
    from redox.core.common import ROOT
    from redox.qm import native_thermal
    original = native_thermal._load_native_surface_record
    body = inspect.getsource(original)
    old = "if r.get('implementation_controls')!=implementation_controls():return None"
    new = "if not _controls_equal(r.get('implementation_controls'),implementation_controls(),ROOT,_prior_root,_receipts):return None"
    if body.count(old) != 1:
        raise ValueError('Native reader changed; inspect before applying path relocation')
    receipts = {}
    namespace = dict(vars(native_thermal), _controls_equal=controls_equal,
                     _prior_root=prior_root, _receipts=receipts)
    local = {}
    exec(compile(body.replace(old, new, 1), __file__, 'exec'), namespace, local)
    native_thermal._load_native_surface_record = local['_load_native_surface_record']
    try:
        yield dict(original_source_root=str(prior_root), relocated_source_root=str(ROOT),
                   native_reader_source_sha256=hashlib.sha256(inspect.getsource(original).encode()).hexdigest(),
                   relocation_module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   controls=receipts, physical_checks_changed=False)
    finally:
        native_thermal._load_native_surface_record = original
