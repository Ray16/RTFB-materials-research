"""Execution-only checkpoints for existing molecular SCF/derivative code.

Class-method wrappers preserve normal binding for copied geometry scanners.
Contexts contain no mean-field references, so completed GPU objects are released.
Objects without a checkpoint context execute their original methods directly.
"""
import json
import os
from pathlib import Path
import time

import numpy as np


def cpu(value):
    return np.asarray(value.get() if hasattr(value,'get') else value)


def write(path,data):
    path=Path(path);temporary=path.with_name('.'+path.name+'.'+str(os.getpid())+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def array(path,value):
    path=Path(path);temporary=path.with_name('.'+path.name+'.'+str(os.getpid())+'.tmp')
    with temporary.open('wb') as stream:np.save(stream,cpu(value),allow_pickle=False)
    temporary.replace(path)


class Context:
    def __init__(self,root,owner):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=False)
        self.owner=owner;self.counter=0;self.current=None

    def begin(self,obj):
        path=self.root/f'{self.counter:04d}';self.counter+=1;path.mkdir();self.current=path
        obj.chkfile=str(path/'SCF.h5')
        write(path/'settings.json',dict(owner=self.owner,coordinates_Bohr=obj.mol.atom_coords().tolist(),
            elements=[obj.mol.atom_symbol(i) for i in range(obj.mol.natm)],charge=obj.mol.charge,
            mult=obj.mol.spin+1,basis=obj.mol.basis,xc=getattr(obj,'xc',None),disp=getattr(obj,'disp',None),
            conv_tol=obj.conv_tol,conv_tol_grad=getattr(obj,'conv_tol_grad',None),max_cycle=obj.max_cycle,
            grid_level=getattr(getattr(obj,'grids',None),'level',None),
            solvent_class=type(obj.with_solvent).__name__ if getattr(obj,'with_solvent',None) else None,
            scientific_settings_changed=False,started=time.time()))
        return path


def wrap_derivative(derivative,filename):
    cls=type(derivative);original=cls.kernel
    if getattr(original,'_redox_checkpoint_wrapper',False):return derivative
    def kernel(self,*args,**kwargs):
        result=original(self,*args,**kwargs)
        context=getattr(self.base,'_redox_checkpoint_context',None)
        if context is not None:
            if context.current is None:raise RuntimeError('Derivative returned before an observed SCF')
            array(context.current/filename,result)
        return result
    kernel._redox_checkpoint_wrapper=True;cls.kernel=kernel
    return derivative


def attach(mf,directory,owner):
    """Record actual solves without binding a scanner to an obsolete MF object."""
    mf._redox_checkpoint_context=Context(directory,owner)
    mf._keys=set(mf._keys)|{'_redox_checkpoint_context'}
    cls=type(mf);original=cls.kernel
    if not getattr(original,'_redox_checkpoint_wrapper',False):
        def kernel(self,*args,**kwargs):
            context=getattr(self,'_redox_checkpoint_context',None)
            if context is None:return original(self,*args,**kwargs)
            path=context.begin(self)
            result=original(self,*args,**kwargs)
            array(path/'density.npy',self.make_rdm1())
            write(path/'SCF_result.json',dict(converged=bool(self.converged),energy_Ha=float(self.e_tot),
                checkpoint=str(path/'SCF.h5'),finished=time.time()))
            return result
        kernel._redox_checkpoint_wrapper=True;cls.kernel=kernel
    for name,filename in [('nuc_grad_method','gradient.npy'),('Hessian','hessian_tensor.npy')]:
        factory=getattr(cls,name)
        if getattr(factory,'_redox_checkpoint_wrapper',False):continue
        def make(self,*args,_factory=factory,_filename=filename,**kwargs):
            result=_factory(self,*args,**kwargs)
            return wrap_derivative(result,_filename) if getattr(self,'_redox_checkpoint_context',None) is not None else result
        make._redox_checkpoint_wrapper=True;setattr(cls,name,make)
    return mf


def recording_builder(builder,root,owner):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    counter=0
    def build(*args,**kwargs):
        nonlocal counter
        directory=root/f'{counter:04d}';counter+=1
        return attach(builder(*args,**kwargs),directory,owner)
    return build


def recording_optimizer(kernel,root,owner):
    """Save optimizer-returned coordinates even when convergence is false."""
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    counter=0
    def run(*args,**kwargs):
        nonlocal counter
        result=kernel(*args,**kwargs);converged,molecule=result
        path=root/f'{counter:04d}.json';counter+=1
        write(path,dict(owner=owner,converged=bool(converged),coordinates_A=molecule.atom_coords(unit='Angstrom').tolist(),
            elements=[molecule.atom_symbol(i) for i in range(molecule.natm)],charge=molecule.charge,
            mult=molecule.spin+1,returned=time.time(),qualified=False))
        return result
    return run
