#!/usr/bin/env python
"""CPU reference honoring the explicitly supplied fixed-complex RI molecule."""
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def molecule_identity(mol):
    """Identify physical/ghost centers and AO basis, independent of range parameter."""
    return (tuple(mol.atom_symbol(i) for i in range(mol.natm)),
            tuple(int(v) for v in mol.atom_charges()),
            mol.atom_coords().tobytes(),repr(mol.basis),int(mol.charge),int(mol.spin),
            int(mol.nelectron),int(mol.nao_nr()))


@contextmanager
def fixed_complex_auxiliary_factory(addons,complex_mol,allowed_molecules):
    """Scope PySCF's rebuilding factory to a fixed complex for the declared states.

    PySCF 2.14 DF.build unconditionally calls make_auxmol, including for the
    range-separated copy. Return a fresh fixed-center auxiliary molecule with
    the caller's range parameter; preserve default behavior for other molecules.
    This process-local context is single-threaded at the Python level and is
    restored even when an SCF fails.
    """
    original=addons.make_auxmol
    allowed={molecule_identity(mol) for mol in allowed_molecules}
    if molecule_identity(complex_mol) not in allowed:raise ValueError('Complex must be in allowed states')
    full_aux=original(complex_mol,None)
    def make_auxmol(mol,auxbasis=None):
        if molecule_identity(mol) not in allowed:return original(mol,auxbasis)
        if auxbasis is not None:raise ValueError('Unexpected non-default auxiliary basis')
        aux=full_aux.copy()
        aux.omega=mol.omega
        return aux
    addons.make_auxmol=make_auxmol
    try:yield
    finally:addons.make_auxmol=original


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(2,minimum_ram_gib=8);os.sched_setaffinity(0,cpu_allocation(2))
    from pyscf.df import addons
    import pyscf
    module_path=ROOT/'scripts/validation/local_counterpoise.py'
    core=api.load_module('counterpoise_fixed_RI_CPU_reference',module_path)
    from pyscf.df.df import DF
    source={str(p):api.digest(p) for p in (Path(__file__).resolve(),module_path,
        Path(__import__('inspect').getfile(DF)),Path(__import__('inspect').getfile(addons)))}
    out=a.out.resolve();provenance=out.with_name(out.name+'_provenance.json')
    if provenance.exists() or out.exists():raise ValueError('Preserve an existing CPU reference')
    record=dict(status='fixed_RI_CPU_reference_running',started=time.time(),source_sha256=source,pyscf_version=pyscf.__version__,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        CPU_affinity=sorted(os.sched_getaffinity(0)),installed_files_modified=False)
    atomic_write_json(provenance,record)
    try:
        coords=[[0,0,0],[0,.76,.58],[0,-.76,.58],[0,0,2.9],[0,.76,3.48],[0,-.76,3.48]]
        molecules=core.molecules(['O','H','H']*2,coords,[0,1,2])
        with fixed_complex_auxiliary_factory(addons,molecules['AB'],molecules.values()):
            result=core.compute(['O','H','H']*2,coords,[0,1,2],'cpu',out)
        api.verify_sources(source);record.update(status='fixed_RI_CPU_reference_computed_review_pending',sums=result['sums'])
    except BaseException as error:record.update(status='fixed_RI_CPU_reference_failed_preserved',error=repr(error));raise
    finally:record['finished']=time.time();atomic_write_json(provenance,record)
    print(json.dumps(result['sums']))


if __name__=='__main__':main()
