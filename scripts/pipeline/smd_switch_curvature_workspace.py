"""Compose bounded response with tiled host switch curvature and stage checkpoints.

No installed file is changed. Full spin/Hessian controls and a separate large
native resource/physical pilot are required before accepting any free energy.
"""
from contextlib import contextmanager
import hashlib
import inspect
from pathlib import Path
import textwrap
import time


def tiled_primitive_source(body):
    old="        p0,p1 = surface['gslice_by_atom'][i_grid_atom]\n"
    if body.count(old)!=1:
        raise ValueError('Installed switch-curvature primitive template changed')
    return body.replace(old,old+'        if p0 == p1:\n            continue\n')


@contextmanager
def compact_smd_execution(directory,budget_GiB):
    import cupy
    import numpy as np
    import psutil
    from gpu4pyscf.solvent.hessian import pcm
    from gpu4pyscf.hessian import rhf
    from gpu4pyscf.solvent import smd
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api

    parent_path=ROOT/'scripts/pipeline/smd_first_response_workspace_v2.py'
    blocks_path=ROOT/'scripts/pipeline/pcm_switch_curvature_blocks.py'
    parent=api.load_module('switch_curvature_response_underlay',parent_path)
    blocks=api.load_module('host_switch_curvature_blocks',blocks_path)
    original=pcm.get_d2F_d2A
    before=inspect.getsource(original);after=tiled_primitive_source(before)
    namespace={};exec(compile(textwrap.dedent(after),str(Path(directory)/'tiled_switch_primitive.py'),'exec'),original.__globals__,namespace)
    primitive=namespace['get_d2F_d2A']
    statistics=dict(status='switch_curvature_execution_active_pending_controls',started=time.time(),
        factories=0,primitive_calls=0,maximum_tile_tensor_bytes=0,host_tensor_bytes=0,
        uploads=0,maximum_uploaded_bytes=0,solute_checkpoints=0,live_binding_refreshes=0,
        grid_tile_points=32,equations_changed=False,surface_grid_changed=False,precision_changed=False)

    def hosted(surface,surface_discretization_method='SWIG'):
        available=psutil.virtual_memory().available
        required=2*len(surface['grid_coords'])*len(surface['atom_coords'])**2*9*8
        if 1.25*required+32*2**30 > .75*available:
            raise MemoryError('Switch-curvature host storage plus headroom is unavailable')
        statistics['factories']+=1
        return blocks.assemble(surface,
            lambda local:primitive(local,surface_discretization_method),cupy,statistics,
            maximum_host_bytes=int(.75*available-32*2**30))

    pcm.get_d2F_d2A=hosted
    try:
        # Patch before the parent captures its transformed-function namespace.
        with parent.compact_smd_execution(directory,budget_GiB=budget_GiB) as proof:
            directory=Path(directory)
            header='\n'.join(Path(pcm.__file__).read_text().splitlines()[:13])+'\n\n'
            (directory/'tiled_switch_primitive.py').write_text(header+after)
            names=('get_dqsym_dx_fix_vgrids','get_dqsym_dx_fix_K_R','analytical_hess_solver')
            previous={name:getattr(pcm,name) for name in names}
            for name,function in previous.items():
                clone=getattr(function,'__kwdefaults__',{}).get('_new')
                if clone is None:
                    raise ValueError('Controlled response wrapper does not expose its cloned function')
                def live(obj,*args,_function=function,_clone=clone,**kwargs):
                    # Native's certified host-curvature context starts later.
                    # Refresh bindings at call time rather than freezing dense
                    # GPU primitives from before that context was installed.
                    for key in ('get_d2F_d2A','get_d2Sii','get_d2D_d2S','get_v_dot_d2A_dot_q'):
                        _clone.__globals__[key]=getattr(pcm,key)
                    statistics['live_binding_refreshes']+=1
                    return _function(obj,*args,**kwargs)
                setattr(pcm,name,live)
            old_kernel=rhf.HessianBase.kernel
            def checkpoint(hess,*args,**kwargs):
                value=old_kernel(hess,*args,**kwargs)
                if isinstance(getattr(hess.base,'with_solvent',None),smd.SMD):
                    index=statistics['solute_checkpoints'];target=directory/f'solute_component_{index:03d}'
                    matrix=value.get() if hasattr(value,'get') else np.asarray(value)
                    np.save(target.with_suffix('.npy'),matrix)
                    atomic_write_json(target.with_suffix('.json'),dict(
                        phase='solute Hessian with equilibrium solvent response; explicit solvent/CDS contributions still pending',
                        coordinates_Bohr=hess.mol.atom_coords().tolist(),basis=repr(hess.mol.basis),
                        charge=hess.mol.charge,mult=hess.mol.spin+1,xc=hess.base.xc,
                        SCF_E_Ha=float(hess.base.e_tot),full_Hessian_qualified=False,G_computed=False))
                    statistics['solute_checkpoints']+=1
                return value
            rhf.HessianBase.kernel=checkpoint
            sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
            statistics['source_sha256']={str(p):sha(p) for p in
                (Path(__file__).resolve(),parent_path,blocks_path,Path(pcm.__file__),Path(rhf.__file__))}
            proof['source_sha256'].update(statistics['source_sha256'])
            proof.update(host_switch_curvature=True,switch_curvature_equations_changed=False)
            try:
                yield proof
                statistics['status']='switch_curvature_finished_requires_controls_and_native_gates'
            except BaseException as error:
                statistics.update(status='switch_curvature_failed_preserved',error=repr(error));raise
            finally:
                rhf.HessianBase.kernel=old_kernel
                for name,function in previous.items():setattr(pcm,name,function)
                statistics['finished']=time.time();atomic_write_json(directory/'switch_curvature.json',statistics)
    finally:
        pcm.get_d2F_d2A=original
