"""Exact ionic composition, spin distinction and actual seeded connectivity."""
import json
from pathlib import Path
from collections import Counter
import numpy as np
from ase import Atoms
from ase.io import read,write
from rdkit import Chem
from pyscf import gto

from redox.core.common import ROOT
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.counterion_geometry import seeds


def test_chloride_microstates_preserve_the_covalent_graft_and_correct_electron_parity(tmp_path):
    package=ROOT/'examples/exact_viologen_corrected_replay/data'
    manifest=json.loads((package/'manifest.json').read_text());parent=next(s for s in manifest['states'] if s['id']=='graft_ox2')
    atoms=read(package/parent['files']['geometry']);smiles=parent['smiles']+'.[Cl-]'
    mol=Chem.MolFromSmiles(smiles);assert Chem.GetFormalCharge(mol)==1
    assert sum(a.GetFormalCharge() for a in mol.GetAtoms() if a.GetAtomicNum()!=17)==2
    cases=list(seeds(atoms));assert len(cases)==7 and len({r[0] for r in cases})==7
    for key,position,definition in cases:
        paired=atoms.copy();paired+=Atoms('Cl',positions=[position]);path=tmp_path/(key+'.xyz');write(path,paired)
        assert hydrogen_topology_status(smiles,path)=='ok'
        assert Counter(paired.get_chemical_symbols())==Counter(atoms.get_chemical_symbols())+Counter({'Cl':1})
        assert min(np.linalg.norm(atoms.positions-position,axis=1))>=2.
        actual=gto.M(atom=list(zip(paired.get_chemical_symbols(),paired.positions.tolist())),basis='def2-svp',charge=1,spin=0,verbose=0)
        assert actual.nelectron%2==0 and actual.nao_nr()==412
    reduced=next(s for s in manifest['states'] if s['id']=='graft_ox1')
    assert reduced['charge']==1 and reduced['mult']==2
