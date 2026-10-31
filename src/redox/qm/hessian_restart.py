"""Frame-checked reuse of completed electronic Hessian components.

Partial arrays remain unqualified until the resumed full calculation passes its
normal derivative/minimum checks. No thermal or electronic energy is borrowed.
"""
from contextlib import contextmanager
import hashlib
import io
import json
from pathlib import Path
from types import MethodType


def response_ceiling_gib(nao,nocc,aligned_atoms=16):
    """Minimum whole-GiB cap for the installed RHF aligned batching estimate."""
    import math
    if any(type(v) is not int or v<1 for v in (nao,nocc,aligned_atoms)) or nocc>nao:
        raise ValueError('Explicit positive AO/occupied/alignment dimensions required')
    minimum=2*aligned_atoms*max(8*3*(nao*nocc*3+nao*nao),8*nao*nao*3*3)
    ceiling=max(3,math.ceil(minimum/2**30))
    if ceiling>24:raise ValueError('Aligned response exceeds this validated 24-GiB execution envelope')
    return ceiling,minimum


def array_digest(array):
    return hashlib.sha256(memoryview(array).cast('B')).hexdigest()


def contiguous_payload(handle,name,content):
    import h5py
    import numpy as np
    dataset=handle[name];creation=dataset.id.get_create_plist();offset=dataset.id.get_offset()
    if (dataset.dtype!=np.dtype('float64') or creation.get_layout()!=h5py.h5d.CONTIGUOUS
            or creation.get_nfilters() or creation.get_external_count() or offset is None or offset<0):
        raise ValueError('Allocated native FP64 contiguous SCF data required')
    count=int(np.prod(dataset.shape));size=count*8
    if offset+size>len(content):raise ValueError('SCF dataset extends beyond its file')
    return np.frombuffer(content,dtype=np.float64,count=count,offset=offset).reshape(dataset.shape)


def read_frame(case):
    """Read completed solute tensors; the outer PCM H1 must be recomputed."""
    import h5py
    import numpy as np
    case=Path(case)
    paths=dict(partial_identity=case/'execution/rks_partial_hess_elec_state_identity.json',
        h1_identity=case/'execution/rks_make_h1_state_identity.json',
        evaluation=case/'points/0000/evaluation.json',settings=case/'SCF_evaluations/0000/0000/settings.json',
        checkpoint=case/'SCF_evaluations/0000/0000/SCF.h5',grids=case/'grids.npz')
    data={k:json.loads(p.read_text()) for k,p in paths.items() if p.suffix=='.json'}
    a,b=data['partial_identity'],data['h1_identity']
    fields=('coordinates_Bohr','basis','charge','mult','xc','SCF_E_Ha','SCF_converged','mo_coeff','mo_energy','mo_occ')
    if any(a[k]!=b[k] for k in fields) or not b['SCF_converged']:
        raise ValueError('Partial Hessian and H1 have different orbital/geometry frames')
    evaluation=data['evaluation'];settings=data['settings']
    coordinates=np.asarray(b['coordinates_Bohr'],dtype=np.float64)
    if (not evaluation['SCF_converged'] or evaluation['coordinates_Bohr']!=b['coordinates_Bohr']
            or abs(evaluation['E_Ha']-b['SCF_E_Ha'])>1e-12
            or settings['coordinates_Bohr']!=b['coordinates_Bohr']
            or settings['charge']!=b['charge'] or settings['mult']!=b['mult']
            or settings['xc']!=b['xc']):
        raise ValueError('Saved SCF and derivative coordinates/state differ')
    arrays={}
    content=paths['checkpoint'].read_bytes()
    with h5py.File(io.BytesIO(content),'r') as handle:
        energy=float(contiguous_payload(handle,'scf/e_tot',content))
        if abs(energy-b['SCF_E_Ha'])>1e-12:raise ValueError('Checkpoint energy differs from the derivative frame')
        for name in ('mo_coeff','mo_occ','mo_energy'):
            value=contiguous_payload(handle,'scf/'+name,content)
            if value.dtype!=np.float64 or list(value.shape)!=b[name]['shape'] or array_digest(value)!=b[name]['sha256']:
                raise ValueError('Checkpoint orbitals do not exactly match the derivative frame')
            arrays[name]=value.copy()
    for name,identity in [('partial',a),('h1',b)]:
        if len(identity['arrays'])!=1:raise ValueError('One complete tensor per derivative stage required')
        path,value=next(iter(identity['arrays'].items()));path=Path(path)
        if not path.is_relative_to(case) or hashlib.sha256(path.read_bytes()).hexdigest()!=value:
            raise ValueError('Changed or foreign derivative tensor')
        paths[name]=path;arrays[name]=np.load(path,mmap_mode='r',allow_pickle=False)
    n=len(coordinates);nao=arrays['mo_coeff'].shape[0];nocc=int((arrays['mo_occ']>0).sum())
    if arrays['partial'].shape!=(n,n,3,3) or arrays['h1'].shape!=(n,3,nao,nocc):
        raise ValueError('Saved component dimensions do not match the declared frame')
    with np.load(paths['grids'],allow_pickle=False) as grids:
        arrays['grid_coords']=grids['coords'];arrays['grid_weights']=grids['weights']
    if (arrays['grid_coords'].shape!=(len(arrays['grid_weights']),3)
            or any(value.dtype!=np.float64 or not np.isfinite(value).all() for value in arrays.values())):
        raise ValueError('Finite FP64 orbital, grid and derivative arrays required')
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()}
    return dict(identity={k:b[k] for k in fields},settings=settings,source_sha256=hashes,
                nocc=nocc,nao=nao,natom=n,H1_scope='solute contribution at the solvated SCF density',
                PCM_H1_included=False),arrays


def restore_mf(mf,frame,arrays,report_path):
    """Restore the original grids/orbitals and verify the resulting energy."""
    import cupy as cp
    from redox.core.io import atomic_write_json
    expected=frame['identity']
    if (mf.mol.atom_coords().tolist()!=expected['coordinates_Bohr'] or mf.mol.charge!=expected['charge']
            or mf.mol.spin+1!=expected['mult'] or mf.xc!=expected['xc']
            or repr(mf.mol.basis)!=expected['basis'] or mf.mol.nao!=frame['nao']):
        raise ValueError('Restart mean-field state or geometry differs')
    for name in ('mo_coeff','mo_occ','mo_energy'):setattr(mf,name,cp.asarray(arrays[name]))
    mf.e_tot=expected['SCF_E_Ha'];mf.converged=True
    mf.grids.coords=cp.asarray(arrays['grid_coords']);mf.grids.weights=cp.asarray(arrays['grid_weights'])
    mf.cphf_grids=mf.grids.copy()
    density=mf.make_rdm1();potential=mf.get_veff(mf.mol,density)
    energy=float(mf.energy_tot(density,mf.get_hcore(),potential));difference=energy-mf.e_tot
    atomic_write_json(report_path,dict(source_energy_Ha=mf.e_tot,recreated_energy_Ha=energy,
        difference_Ha=difference,original_grid_arrays_available=True,native_G_qualified=False))
    if abs(difference)>1e-8:raise ValueError('Exact-grid restarted density energy does not reproduce')


@contextmanager
def reuse_components(hessian,arrays,directory):
    """Reuse solute stages while preserving the outer PCM/SMD H1 addition."""
    import cupy as cp
    from redox.core.io import atomic_write_json
    directory=Path(directory);directory.mkdir(exist_ok=False)
    if not hasattr(hessian,'undo_solvent') or len(hessian.__class__.__bases__)!=2:
        raise ValueError('Explicit two-base solvent Hessian wrapper required')
    base_class=hessian.__class__.__bases__[1]
    old_partial=('partial_hess_elec' in hessian.__dict__,hessian.partial_hess_elec)
    old_h1=('make_h1' in base_class.__dict__,base_class.make_h1)
    counts=dict(partial_hess_elec=0,make_h1=0)
    try:
        def partial(self,*args,**kwargs):
            counts['partial_hess_elec']+=1
            return cp.asarray(arrays['partial']).copy()
        def solute_h1(self,*args,**kwargs):
            if self is not hessian:return old_h1[1](self,*args,**kwargs)
            counts['make_h1']+=1
            # The outer WithSolventHess.make_h1 still adds its PCM term to
            # this fresh copy. The stored solute checkpoint stays read-only.
            return cp.asarray(arrays['h1']).copy()
        hessian.partial_hess_elec=MethodType(partial,hessian)
        base_class.make_h1=solute_h1
        yield counts
        if any(v!=1 for v in counts.values()):raise ValueError('Each exact saved stage must be consumed once')
    finally:
        if old_partial[0]:hessian.partial_hess_elec=old_partial[1]
        else:delattr(hessian,'partial_hess_elec')
        if old_h1[0]:base_class.make_h1=old_h1[1]
        else:delattr(base_class,'make_h1')
        atomic_write_json(directory/'reuse.json',dict(stage_calls=counts,checkpoint_arrays_modified=False,
            H1_checkpoint_scope='solute contribution at the solvated SCF density',PCM_H1_recomputed=True,
            original_grid_arrays_available=True,full_physical_QC_required=True,native_G_qualified=False))


@contextmanager
def capture_solution_h1(hessian,directory,identity):
    """Save the complete solvent-wrapped H1 after its actual computation."""
    import os
    import numpy as np
    from redox.core.io import atomic_write_json
    directory=Path(directory);directory.mkdir(exist_ok=False)
    present='make_h1' in hessian.__dict__;original=hessian.make_h1;calls=0
    def capture(self,*args,**kwargs):
        nonlocal calls
        value=original(*args,**kwargs);calls+=1
        if calls!=1:raise ValueError('One complete RKS solution H1 per control required')
        array=np.asarray(value.get() if hasattr(value,'get') else value)
        if array.dtype!=np.float64 or not np.isfinite(array).all():raise ValueError('Finite FP64 full solution H1 required')
        path=directory/'full_H1.npy';temporary=directory/'full_H1.npy.partial'
        with temporary.open('xb') as stream:np.save(stream,array,allow_pickle=False);stream.flush();os.fsync(stream.fileno())
        os.replace(temporary,path)
        atomic_write_json(directory/'identity.json',dict(scope='full SMD H1 including PCM contribution',identity=identity,
            shape=list(array.shape),dtype='float64',path=str(path),file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            array_sha256=array_digest(array),native_G_qualified=False))
        return value
    hessian.make_h1=MethodType(capture,hessian)
    try:yield
    finally:
        if present:hessian.make_h1=original
        else:delattr(hessian,'make_h1')
