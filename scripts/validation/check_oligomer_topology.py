#!/usr/bin/env python
"""Check graph topology and its realization in every available model XYZ."""
import argparse
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.geometry import _skeleton
from redox.build.oligomer_topology import chain_topology


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    feed=json.loads((ROOT/'results/polymer_models/interaction_checked.json').read_text());results=[]
    for row in feed['models']:
        declared=Chem.MolFromSmiles(row['smiles']);parts=Chem.GetMolFrags(declared,asMols=True)
        expected_fragments=2 if row['model_type']=='two_nonbonded_segments' else 1
        assert len(parts)==expected_fragments
        topology=[chain_topology(Chem.MolToSmiles(part,isomericSmiles=True),1 if expected_fragments==2 else row['repeats']) for part in parts]
        paths=[('starting',row['input_xyz'])]+([('completed_DFT',row['final_xyz'])] if row['final_xyz'] else [])
        for stage,name in paths:
            path=Path(name);observed=Chem.MolFromXYZFile(str(path));assert observed is not None
            rdDetermineBonds.DetermineConnectivity(observed,charge=0);reference=Chem.AddHs(declared)
            mapping=_skeleton(observed).GetSubstructMatch(_skeleton(reference))
            assert len(mapping)==reference.GetNumAtoms() and observed.GetNumAtoms()==reference.GetNumAtoms(),row['key']
            assert _skeleton(reference).HasSubstructMatch(_skeleton(observed)),row['key']
            results.append(dict(key=row['key'],stage=stage,graph_topology=topology,actual_XYZ_matches_full_graph=True,
                                xyz=str(path),xyz_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    atomic_write_json(args.out,dict(status='checked',models=len(feed['models']),starting_XYZ_checked=sum(r['stage']=='starting' for r in results),
        completed_DFT_XYZ_checked=sum(r['stage']=='completed_DFT' for r in results),
        scope='Open 2n-carbon chain, alternating phenyl sites, para benzyl linkers, fragment count and complete XYZ connectivity; not a free-energy/minimum test',results=results))
    print('Topology verified for',len(results),'model geometries')


if __name__=='__main__':main()
