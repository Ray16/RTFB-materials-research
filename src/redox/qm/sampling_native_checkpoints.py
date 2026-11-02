"""Full-orbital native SCF evidence and center grids without extra response cubes.

The native solver keeps all virtual orbitals. Only the already-written partial
Hessian/H1 and exact center SCF/grid frame are restart inputs; no complete native
G or response restart is inferred from these files alone.
"""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path

import numpy as np


def cpu(value):
    return np.asarray(value.get() if hasattr(value,'get') else value)


def write(path,value):
    temporary=path.with_name('.'+path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def capture_before_hessian(mf, root):
    root=Path(root);context=mf._redox_checkpoint_context
    if context.current is None or not mf.converged:
        raise ValueError('Converged actual center checkpoint required before a native Hessian')
    checkpoint=context.current/'SCF.h5'
    if not checkpoint.is_file():raise ValueError('Durable actual center SCF is missing')
    coords=cpu(mf.grids.coords);weights=cpu(mf.grids.weights)
    matched=getattr(mf,'cphf_grids',None)
    if (matched is None or not np.array_equal(cpu(matched.coords),coords)
            or not np.array_equal(cpu(matched.weights),weights)
            or coords.shape!=(len(weights),3) or not np.isfinite(coords).all() or not np.isfinite(weights).all()):
        raise ValueError('Exact matched finite SCF/CPHF grids required')
    with (root/'.grids.tmp').open('wb') as stream:np.savez(stream,coords=coords,weights=weights)
    (root/'.grids.tmp').replace(root/'grids.npz')
    write(root/'restart_frame.json',dict(layout='sampling_native.v1',SCF_checkpoint=str(checkpoint),
        checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        geometry=str(root/'geometry.xyz'),coordinates_Bohr=mf.mol.atom_coords().tolist(),
        charge=mf.mol.charge,mult=mf.mol.spin+1,basis=mf.mol.basis,xc=mf.xc,
        SCF_E_Ha=float(mf.e_tot),SCF_converged=True,
        grids=str(root/'grids.npz'),grids_sha256=hashlib.sha256((root/'grids.npz').read_bytes()).hexdigest(),
        CPHF_grids_matched=True,partial_tensors_location=str(root/'execution'),
        raw_response_vectors_retained=False,native_G_qualified=False))


@contextmanager
def recording_native(dft, checkpoints, root, owner):
    original=dft._build_smd_mf;root=Path(root);counter=0
    def build(*args,**kwargs):
        nonlocal counter
        folder=root/'SCF_frames'/f'{counter:04d}';counter+=1
        mf=checkpoints.attach(original(*args,**kwargs),folder,owner,save_density=False)
        mf._redox_sampling_native_root=str(root);mf._keys=set(mf._keys)|{'_redox_sampling_native_root'}
        cls=type(mf);factory=cls.Hessian
        if not getattr(factory,'_redox_native_grid_wrapper',False):
            def make(self,*args,**kwargs):
                output=getattr(self,'_redox_sampling_native_root',None)
                if output is not None:capture_before_hessian(self,output)
                return factory(self,*args,**kwargs)
            make._redox_native_grid_wrapper=True;make._redox_checkpoint_wrapper=True;cls.Hessian=make
        return mf
    dft._build_smd_mf=build
    try:yield
    finally:dft._build_smd_mf=original
