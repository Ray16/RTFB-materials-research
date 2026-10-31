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


def inspect_source(case,frame):
    """Bind raw solver blocks to an explicitly complete solvent H1/frame.

    Closure/liveness is a workflow admission check. This function reads only
    immutable completed arrays and does not qualify a Hessian or free energy.
    """
    import numpy as np
    case=Path(case).resolve();plan=json.loads((case/'plan.json').read_text())
    if (plan.get('role')!='screening_polymer_Hessian_control' or plan.get('screening_variant')!='corrected'
            or not plan.get('full_polymer_Hessian_control')):
        raise ValueError('Only a full corrected solvent-wrapped Hessian may supply response blocks')
    center=case/'points/0000/evaluation.json';point=json.loads(center.read_text());grids=case/'grids.npz'
    if (not point['SCF_converged'] or point['coordinates_Bohr']!=frame['identity']['coordinates_Bohr']
            or not math.isfinite(point['E_Ha']) or abs(point['E_Ha']-frame['identity']['SCF_E_Ha'])>1e-12):
        raise ValueError('Response center differs from the original orbital frame')
    identity_path=case/'solution_H1/identity.json';capture=json.loads(identity_path.read_text())
    expected=dict(center_evaluation=str(center),center_evaluation_sha256=digest(center),grids_sha256=digest(grids),
                  source_namespace=plan['original_namespace'],screening_variant='corrected')
    path=Path(capture['path']).resolve();array=np.load(path,mmap_mode='r',allow_pickle=False)
    if (capture['scope']!='full SMD H1 including PCM contribution' or capture['identity']!=expected
            or not path.is_relative_to(case/'solution_H1') or digest(path)!=capture['file_sha256']
            or array.dtype!=np.float64 or list(array.shape)!=capture['shape']
            or array.shape!=(frame['natom'],3,frame['nao'],frame['nocc']) or not np.isfinite(array).all()
            or hashlib.sha256(memoryview(array).cast('B')).hexdigest()!=capture['array_sha256']):
        raise ValueError('Complete, source-bound finite FP64 solution H1 required')
    paths=[case/'plan.json',center,grids,identity_path,path];blocks=sorted((case/'response_blocks').glob('block_*.json'))
    if not blocks:raise ValueError('No complete raw response block is available')
    total=0;solver_shas=set()
    for index,metadata in enumerate(blocks):
        if metadata.name!=f'block_{index:03d}.json':raise ValueError('Response blocks must form a contiguous completed prefix')
        row=json.loads(metadata.read_text());saved=row['checkpoint'];payload=Path(saved['path']).resolve()
        value=np.load(payload,mmap_mode='r',allow_pickle=False)
        if (row['kind']!='raw_Krylov_solution_before_occupied_block_reconstruction' or row['identity']!=expected
                or not payload.is_relative_to(case/'response_blocks') or digest(payload)!=saved['file_sha256']
                or value.dtype!=np.float64 or not np.isfinite(value).all() or list(value.shape)!=saved['shape']
                or value.ndim!=2 or value.shape[1]!=frame['nao']*frame['nocc']
                or hashlib.sha256(memoryview(value).cast('B')).hexdigest()!=saved['array_sha256']):
            raise ValueError('Changed or mismatched source response payload')
        total+=value.shape[0];solver_shas.add(row['original_solver_source_sha256']);paths.extend([metadata,payload])
    if total>3*frame['natom'] or len(solver_shas)!=1:raise ValueError('Response prefix dimensions or solver identity differ')
    return dict(source_case=str(case),blocks=str(case/'response_blocks'),identity=expected,full_H1=capture,
        completed_blocks=len(blocks),completed_RHS=total,original_solver_source_sha256=next(iter(solver_shas)),
        source_sha256={str(p):digest(p) for p in paths},native_G_qualified=False)


@contextmanager
def checked_blocks(module,source,directory,*,expected_identity=None,source_h1=None,current_h1=None):
    import numpy as np
    from redox.core.io import atomic_write_json
    source=Path(source);directory=Path(directory);directory.mkdir(exist_ok=False)
    original=module.residual_krylov;original_sha=getattr(original,'_redox_original_solver_source_sha256',digest(original.__code__.co_filename));calls=[]
    def solve(aop,b,*args,**kwargs):
        index=len(calls);entry=dict(index=index,reused=False);calls.append(entry)
        on_gpu=type(b).__module__.split('.')[0]=='cupy'
        if on_gpu:
            import cupy as cp
            host=cp.asnumpy;device=cp.asarray
        else:host=np.asarray;device=np.asarray
        rhs=np.asarray(host(b));metadata=source/f'block_{index:03d}.json'
        if rhs.dtype!=np.float64 or not np.isfinite(rhs).all():raise ValueError('Finite FP64 response RHS required')
        if metadata.exists() and not args:
            old=json.loads(metadata.read_text());checkpoint=old['checkpoint'];path=Path(checkpoint['path'])
            if expected_identity is not None and old.get('identity')!=expected_identity:raise ValueError('Saved response metadata changed its bound frame')
            if (not path.resolve().is_relative_to(source.resolve()) or digest(path)!=checkpoint['file_sha256']):
                raise ValueError('Changed or foreign saved response block')
            candidate=np.load(path,allow_pickle=False)
            if (candidate.dtype!=np.float64 or not np.isfinite(candidate).all()
                    or list(candidate.shape)!=checkpoint['shape']
                    or hashlib.sha256(candidate.tobytes()).hexdigest()!=checkpoint['array_sha256']):
                raise ValueError('Invalid saved FP64 response payload')
            same=(old['rhs_sha256']==hashlib.sha256(rhs.tobytes()).hexdigest()
                  and old['original_solver_source_sha256']==original_sha and candidate.shape==rhs.shape)
            if source_h1 is not None:
                actual=json.loads(Path(current_h1).read_text())
                if actual.get('equilibrium_solvation') is not True:raise ValueError('Fresh equilibrium-solvent context must be witnessed before response reuse')
                same=(same and actual['scope']==source_h1['scope']=='full SMD H1 including PCM contribution'
                      and actual['array_sha256']==source_h1['array_sha256']
                      and actual['identity']['grids_sha256']==source_h1['identity']['grids_sha256']
                      and actual['identity']['source_namespace']==source_h1['identity']['source_namespace'])
            if same:
                # Some response implementations reuse their input as work
                # space. Keep the candidate returned to the caller unchanged.
                value=device(candidate);response=np.asarray(host(aop(value.copy())))
                if response.shape!=rhs.shape or response.dtype!=np.float64 or not np.isfinite(response).all():raise ValueError('Fresh FP64 response operator failed')
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
            else:entry['reason']='RHS_shape_full_H1_or_solver_identity_differs_recompute'
        result=original(aop,b,*args,**kwargs)
        atomic_write_json(directory/'status.json',dict(calls=calls,native_G_qualified=False))
        return result
    solve._redox_original_solver_source_sha256=original_sha
    module.residual_krylov=solve
    try:yield calls
    finally:
        module.residual_krylov=original
        atomic_write_json(directory/'status.json',dict(calls=calls,source=str(source),
            original_solver_sha256=original_sha,equations_changed=False,native_G_qualified=False))
