#!/usr/bin/env python
"""Add durable SCF files to an owned Python worker without changing its solver."""
import argparse
from contextlib import contextmanager
import os
from pathlib import Path
import runpy
import sys
import time

from redox.core.io import atomic_write_json


@contextmanager
def checkpoints(directory):
    from redox.qm import dft
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=False)
    original=dft._build_smd_mf
    count=0
    def build(*args,**kwargs):
        nonlocal count
        mf=original(*args,**kwargs)
        folder=directory/f'SCF_{count:04d}'
        count+=1
        folder.mkdir()
        mf.chkfile=str(folder/'SCF.h5')
        settings=dict(charge=mf.mol.charge,mult=mf.mol.spin+1,basis=repr(mf.mol.basis),
            xc=mf.xc,disp=getattr(mf,'disp',None),coordinates_Bohr=mf.mol.atom_coords().tolist(),
            solvent=mf.with_solvent.solvent,started=time.time(),qualified_result=False)
        atomic_write_json(folder/'settings.json',settings)
        previous=getattr(mf,'callback',None)
        def callback(environment):
            if callable(previous):previous(environment)
            atomic_write_json(folder/'progress.json',dict(cycle=environment.get('cycle'),
                E_Ha=float(environment['e_tot']),updated=time.time(),qualified_result=False))
        mf.callback=callback
        return mf
    dft._build_smd_mf=build
    try:yield
    finally:dft._build_smd_mf=original


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('program',type=Path)
    parser.add_argument('arguments',nargs=argparse.REMAINDER)
    args=parser.parse_args()
    before=sys.argv
    try:
        with checkpoints(args.out):
            sys.argv=[str(args.program),*args.arguments]
            runpy.run_path(str(args.program),run_name='__main__')
    finally:sys.argv=before


if __name__=='__main__':main()
