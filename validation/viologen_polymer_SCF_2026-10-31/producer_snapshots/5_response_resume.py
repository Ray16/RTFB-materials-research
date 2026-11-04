"""Reuse a saved linear-response candidate only after a fresh equation check.

The caller must bind the same orbital/grid/full-H1/operator frame. This layer
also checks the actual RHS bytes and the original solver source. Invalid or
corrupt backing records are never silently admitted.
"""
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@contextmanager
def checked_blocks(module,source,directory):
    import numpy as np
    from redox.core.io import atomic_write_json
    source=Path(source);directory=Path(directory);directory.mkdir(exist_ok=False)
    original=module.residual_krylov;original_sha=digest(original.__code__.co_filename);calls=[]
    def solve(aop,b,*args,**kwargs):
        index=len(calls);entry=dict(index=index,reused=False);calls.append(entry)
        on_gpu=type(b).__module__.split('.')[0]=='cupy'
        if on_gpu:
            import cupy as cp
            host=cp.asnumpy;device=cp.asarray
        else:host=np.asarray;device=np.asarray
        rhs=np.asarray(host(b));metadata=source/f'block_{index:03d}.json'
        if metadata.exists() and not args:
            old=json.loads(metadata.read_text());checkpoint=old['checkpoint'];path=Path(checkpoint['path'])
            if (not path.is_relative_to(source) or digest(path)!=checkpoint['file_sha256']):
                raise ValueError('Changed or foreign saved response block')
            candidate=np.load(path,allow_pickle=False)
            if (candidate.dtype!=np.float64 or not np.isfinite(candidate).all()
                    or list(candidate.shape)!=checkpoint['shape']
                    or hashlib.sha256(candidate.tobytes()).hexdigest()!=checkpoint['array_sha256']):
                raise ValueError('Invalid saved FP64 response payload')
            same=(old['rhs_sha256']==hashlib.sha256(rhs.tobytes()).hexdigest()
                  and old['original_solver_source_sha256']==original_sha and candidate.shape==rhs.shape)
            if same:
                value=device(candidate);response=np.asarray(host(aop(value)))
                if response.shape!=rhs.shape or not np.isfinite(response).all():raise ValueError('Fresh response operator failed')
                residual=candidate+response-rhs;norm=float(np.linalg.norm(residual))
                tolerance=kwargs.get('tol',1e-10);target=tolerance*.01
                entry.update(source=str(metadata),fresh_true_residual_norm=norm,current_absolute_target=target)
                if math.isfinite(norm) and norm<=target:
                    record=kwargs.get('record')
                    if not callable(record):raise ValueError('Original source-bound response recorder required')
                    record(dict(algorithm='Saved candidate, freshly checked unchanged I+A equation',
                        requested_absolute_tolerance=tolerance,solver_absolute_target=target,true_residual_norm=norm,
                        maximum_absolute_residual=float(abs(residual).max()),rhs_norm=float(np.linalg.norm(rhs)),
                        solver_info=0,operator_evaluations=1,iterations=0,residual_history=[],shape=list(rhs.shape),
                        precision='float64',level_shift_added=False,checkpoint_reused=True))
                    entry['reused']=True;atomic_write_json(directory/'status.json',dict(calls=calls,native_G_qualified=False))
                    return value
                entry['reason']='candidate_failed_fresh_residual_check_recompute_original_equation'
            else:entry['reason']='RHS_shape_or_solver_identity_differs_recompute'
        result=original(aop,b,*args,**kwargs)
        atomic_write_json(directory/'status.json',dict(calls=calls,native_G_qualified=False))
        return result
    module.residual_krylov=solve
    try:yield calls
    finally:
        module.residual_krylov=original
        atomic_write_json(directory/'status.json',dict(calls=calls,source=str(source),
            original_solver_sha256=original_sha,equations_changed=False,native_G_qualified=False))
