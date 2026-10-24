"""Stream exact analytic PCM curvature contractions on bounded GPU tiles.

The installed FP64 primitives, every atom and every surface point are retained.
This profile requires full spin/Hessian replay before production use. It changes
contraction/storage order, never a functional, grid, derivative or QC threshold.
"""
from collections import OrderedDict
from contextlib import contextmanager
import inspect
import math
from pathlib import Path
import textwrap
import time

import numpy as np


def vectorized_switch_source(body):
    loop = ('            for i_atom in range(natom):\n'
            '                d2fiJK[:, i_atom, i_atom, :, :] = d2fiJK_diagonal[:, i_atom, :, :]')
    if body.count(loop) != 2:
        raise ValueError('Installed switch-curvature diagonal template changed')
    body = body.replace(loop, '            diagonal = cupy.arange(natom)\n'
        '            d2fiJK[:, diagonal, diagonal, :, :] = d2fiJK_diagonal')
    owner = "        p0,p1 = surface['gslice_by_atom'][i_grid_atom]\n"
    if body.count(owner) != 1:
        raise ValueError('Installed grid ownership template changed')
    return body.replace(owner, owner+'        if p0 == p1:\n            continue\n')


class SwitchSource:
    def __init__(self, surface, primitive, xp, metrics, progress, tile_points=32):
        self.surface, self.primitive, self.xp = surface, primitive, xp
        self.metrics, self.progress, self.tile = metrics, progress, tile_points
        self.points, self.atoms = len(surface['grid_coords']), len(surface['atom_coords'])

    def contract(self, kind, vector):
        from redox.core.common import ROOT
        from redox.workflows import polymer_campaign as api
        blocks = api.load_module('streamed_PCM_surface_slices', ROOT/'scripts/pipeline/pcm_switch_curvature_blocks.py')
        xp = self.xp
        if kind not in ('F','A'):
            raise ValueError('Unknown switch-curvature kind')
        if vector.shape != (self.points,) or vector.dtype != xp.float64:
            raise ValueError('FP64 surface-vector contraction required')
        result = xp.zeros((self.atoms,self.atoms,3,3),dtype=xp.float64)
        started = time.time()
        self.metrics['switch_contractions'] += 1
        index = self.metrics['switch_contractions']
        for begin in range(0,self.points,self.tile):
            end = min(self.points,begin+self.tile)
            values = self.primitive(blocks.local_surface(self.surface,begin,end),'SWIG')
            selected = values[0 if kind == 'F' else 1]
            result += xp.einsum('abxyq,q->abxy',selected,vector[begin:end])
            self.metrics['switch_primitive_calls'] += 1
            self.metrics['maximum_switch_tile_bytes'] = max(self.metrics['maximum_switch_tile_bytes'],sum(v.nbytes for v in values))
            if begin == 0 or end == self.points or (begin//self.tile+1) % 64 == 0:
                self.progress(index,kind,end,self.points,result,started)
            del selected, values
        return result


class WeightedSwitch:
    def __init__(self, source, kind):
        self.source, self.kind = source, kind
        self.shape = (source.atoms,source.atoms,3,3,source.points)
        self.ndim, self.dtype = 5, np.dtype('float64')

    def __matmul__(self, vector):
        return self.source.contract(self.kind,vector)


class StreamedDiagonal:
    def __init__(self,surface,first,second,xp,metrics):
        if not isinstance(second,WeightedSwitch) or second.kind != 'F':
            raise ValueError('Diagonal operator requires its streamed switch parent')
        self.surface,self.first,self.second = surface,first,second
        self.xp,self.metrics = xp,metrics

    def __matmul__(self,vector):
        xp = self.xp
        derivative = xp.ascontiguousarray(self.first.transpose(2,0,1))
        atoms,axes,points = derivative.shape
        if axes != 3 or vector.shape != (points,):
            raise ValueError('Unexpected diagonal PCM contraction shape')
        switch,charge = self.surface['switch_fun'],self.surface['charge_exp']
        coefficient = math.sqrt(2/math.pi)*charge*vector
        flat = derivative.reshape(atoms*3,points)
        first = ((flat*(2*coefficient/switch**3))@flat.T).reshape(atoms,3,atoms,3).transpose(0,2,1,3)
        self.metrics['diagonal_contractions'] += 1
        return first-self.second@(coefficient/switch**2)


class GPUCurvatureBlocks:
    def __init__(self,surface,primitive,xp,metrics,cache_bytes=128*1024**2):
        self.surface,self.primitive,self.xp,self.metrics = surface,primitive,xp,metrics
        self.points = len(surface['charge_exp'])
        self.cache,self.bytes,self.limit = OrderedDict(),0,cache_bytes

    def block(self,rows,columns):
        xp = self.xp
        key = (tuple(rows),tuple(columns))
        if key in self.cache:
            self.cache.move_to_end(key)
            self.metrics['offdiagonal_cache_hits'] += 1
            return self.cache[key]
        if len(rows) == 0 or len(columns) == 0:
            return tuple(xp.zeros((3,3,len(rows),len(columns)),dtype=xp.float64) for _ in range(2))
        combined = np.concatenate([rows,columns])
        indices,inverse = np.unique(combined,return_inverse=True)
        expected_bytes = 2*9*len(indices)**2*8
        if expected_bytes > 512*1024**2:
            raise MemoryError('PCM primitive union exceeds the bounded GPU tile')
        local = {name:self.surface[name][xp.asarray(indices)] for name in ('charge_exp','grid_coords','norm_vec')}
        D,S = self.primitive(local,with_S=True,with_D=True)
        left,right = xp.asarray(inverse[:len(rows)]),xp.asarray(inverse[len(rows):])
        values = tuple(xp.ascontiguousarray(xp.take(xp.take(a,left,axis=2),right,axis=3)) for a in (D,S))
        self.metrics['offdiagonal_primitive_calls'] += 1
        self.metrics['maximum_offdiagonal_primitive_bytes'] = max(self.metrics['maximum_offdiagonal_primitive_bytes'],D.nbytes+S.nbytes)
        size = sum(a.nbytes for a in values)
        if size <= self.limit:
            while self.cache and self.bytes+size > self.limit:
                _,old = self.cache.popitem(last=False)
                self.bytes -= sum(a.nbytes for a in old)
            self.cache[key] = values
            self.bytes += size
            self.metrics['peak_GPU_block_cache_bytes'] = max(self.metrics['peak_GPU_block_cache_bytes'],self.bytes)
        return values


class GPUOffDiagonal:
    def __init__(self,blocks,kind,transposed=False):
        self.blocks,self.kind,self.transposed = blocks,kind,transposed

    def transpose(self,*axes):
        if axes != (0,1,3,2):
            raise ValueError('Unsupported PCM surface-axis transpose')
        return GPUOffDiagonal(self.blocks,self.kind,not self.transposed)

    def __getitem__(self,key):
        if len(key) != 4:
            raise ValueError('Expected Cartesian/Cartesian/surface/surface slice')
        rows,columns = np.arange(self.blocks.points)[key[2]],np.arange(self.blocks.points)[key[3]]
        if np.ndim(rows) != 1 or np.ndim(columns) != 1:
            raise ValueError('Scalar surface indexing is unsupported')
        original_rows,original_columns = (columns,rows) if self.transposed else (rows,columns)
        value = self.blocks.block(original_rows,original_columns)[self.kind]
        if self.transposed:
            value = value.transpose(0,1,3,2)
        self.blocks.metrics['offdiagonal_block_requests'] += 1
        return self.blocks.xp.ascontiguousarray(value[key[0],key[1]])


@contextmanager
def compact_smd_execution(directory,budget_GiB):
    import cupy as cp
    from gpu4pyscf.solvent.hessian import pcm
    from gpu4pyscf.hessian import rhf
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api
    directory = Path(directory)
    underlay = api.load_module('streamed_PCM_first_response',ROOT/'scripts/pipeline/smd_first_response_workspace_v2.py')
    metrics = dict(status='streamed_PCM_active_unqualified',started=time.time(),switch_contractions=0,
        switch_primitive_calls=0,maximum_switch_tile_bytes=0,diagonal_contractions=0,
        offdiagonal_primitive_calls=0,offdiagonal_cache_hits=0,offdiagonal_block_requests=0,
        maximum_offdiagonal_primitive_bytes=0,peak_GPU_block_cache_bytes=0,
        precision='float64',surface_grid_changed=False,atoms_omitted=False,
        full_spin_replay_required=True,large_resource_qualification_required=True)
    with underlay.compact_smd_execution(directory,budget_GiB=budget_GiB) as proof:
        names=('get_d2F_d2A','get_d2Sii','get_d2D_d2S','get_v_dot_d2A_dot_q')
        originals={n:getattr(pcm,n) for n in names}
        clones=[]
        for name in ('get_dqsym_dx_fix_vgrids','get_dqsym_dx_fix_K_R','analytical_hess_solver'):
            clone=getattr(getattr(pcm,name),'__kwdefaults__',{}).get('_new')
            if clone is None:
                raise ValueError('Controlled response must expose its function namespace')
            clones.append(clone)
        body=vectorized_switch_source(inspect.getsource(originals['get_d2F_d2A']))
        scope=dict(originals['get_d2F_d2A'].__globals__)
        exec(compile(textwrap.dedent(body),str(directory/'vectorized_switch_primitive.py'),'exec'),scope)
        primitive=scope['get_d2F_d2A']
        header='\n'.join(Path(pcm.__file__).read_text().splitlines()[:13])+'\n\n'
        (directory/'vectorized_switch_primitive.py').write_text(header+body)

        def progress(index,kind,completed,total,value,started):
            # Partial contractions remain unqualified. Small dense outputs are
            # durable even if a later native stage is interrupted.
            target=directory/f'switch_{index:03d}_{completed:06d}.npy'
            np.save(target,cp.asnumpy(value))
            atomic_write_json(directory/'stream_progress.json',dict(contraction=index,kind=kind,
                completed_surface_points=completed,total_surface_points=total,elapsed_seconds=time.time()-started,
                checkpoint=str(target),checkpoint_sha256=api.digest(target),qualified_native_G=False))

        def switches(surface,surface_discretization_method='SWIG'):
            if surface_discretization_method.upper()!='SWIG':
                raise ValueError('This profile is qualified only for SWIG SMD surfaces')
            source=SwitchSource(surface,primitive,cp,metrics,progress)
            return WeightedSwitch(source,'F'),WeightedSwitch(source,'A')

        def diagonal(surface,first,second,stream=None):
            return StreamedDiagonal(surface,first,second,cp,metrics)

        def offdiagonal(surface,with_S=True,with_D=False,stream=None):
            blocks=GPUCurvatureBlocks(surface,originals['get_d2D_d2S'],cp,metrics)
            return GPUOffDiagonal(blocks,0) if with_D else None,GPUOffDiagonal(blocks,1) if with_S else None

        pcm.get_d2F_d2A,pcm.get_d2Sii,pcm.get_d2D_d2S=switches,diagonal,offdiagonal
        pcm.get_v_dot_d2A_dot_q=lambda tensor,left,right:tensor@(left*right)
        for clone in clones:
            for key in names:
                clone.__globals__[key]=getattr(pcm,key)
        original_kernel=rhf.HessianBase.kernel
        def checkpoint(hess,*args,**kwargs):
            value=original_kernel(hess,*args,**kwargs)
            array=value.get() if hasattr(value,'get') else np.asarray(value)
            np.save(directory/'solute_component.npy',array)
            atomic_write_json(directory/'solute_component.json',dict(coordinates_Bohr=hess.mol.atom_coords().tolist(),
                basis=repr(hess.mol.basis),charge=hess.mol.charge,mult=hess.mol.spin+1,
                xc=hess.base.xc,SCF_E_Ha=float(hess.base.e_tot),qualified_native_G=False))
            return value
        rhf.HessianBase.kernel=checkpoint
        metrics['source_sha256']={str(p):api.digest(p) for p in (Path(__file__).resolve(),Path(pcm.__file__),ROOT/'scripts/pipeline/pcm_switch_curvature_blocks.py')}
        proof.update(streamed_PCM_curvature=True,host_full_curvature_tensors=False,
            maximum_GPU_curvature_cache_bytes=128*1024**2)
        proof['source_sha256'].update(metrics['source_sha256'])
        try:
            yield proof
            metrics['status']='streamed_PCM_finished_requires_full_replay_and_physical_QC'
        except BaseException as error:
            metrics.update(status='streamed_PCM_failed_preserved',error=repr(error))
            raise
        finally:
            rhf.HessianBase.kernel=original_kernel
            for name,value in originals.items():setattr(pcm,name,value)
            metrics['finished']=time.time()
            atomic_write_json(directory/'streamed_PCM.json',metrics)
