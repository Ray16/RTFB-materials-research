#!/usr/bin/env python
"""Audit geometry diversity of declared starts without implying minima or weights."""
import argparse
import json
import os
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def separated_subset(matrix,cutoff):
    matrix=np.asarray(matrix,dtype=float)
    if (matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or not np.isfinite(matrix).all()
            or np.min(matrix)<0 or np.max(abs(matrix-matrix.T))>1e-8):
        raise ValueError('Finite symmetric nonnegative RMSD matrix required')
    chosen=[]
    for index in range(len(matrix)):
        if all(matrix[index,old]>=cutoff for old in chosen):chosen.append(index)
    return chosen


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,action='append',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    from rdkit import Chem,rdBase
    from rdkit.Chem import rdMolAlign
    helper_path=ROOT/'scripts/pipeline/polymer_conformer_thermals.py'
    helper=api.load_module('seed_audit_declared_geometry_mapping',helper_path)
    sources={str(p):api.digest(p) for p in (Path(__file__).resolve(),helper_path)}
    rows=[];args.out.parent.mkdir(parents=True,exist_ok=True)
    for root in args.root:
        for path in sorted(root.resolve().glob('*/plan.json')):
            plan=api.read_plan(path.parent)
            if len(plan.get('tasks',[]))!=20 or not plan['tasks'][0]['key'].startswith(('all_same__','alternating__')):continue
            api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
            tasks=plan['tasks'];smiles={t['smiles'] for t in tasks}
            if len(smiles)!=1:raise ValueError('One exact graph per conformer batch required')
            molecules=[]
            for task in tasks:
                xyz=path.parent/'live_inputs'/task['key']/'input.xyz'
                molecules.append(helper.alignment_molecule(task['smiles'],xyz));sources[str(xyz)]=api.digest(xyz)
            template=molecules[0]
            matches=template.GetSubstructMatches(template,uniquify=False,useChirality=True,maxMatches=10001)
            if len(matches)>=10001:raise ValueError('Symmetry enumeration limit reached; diversity needs separate review')
            maps=[list(enumerate(match)) for match in matches]
            matrix=np.zeros((20,20))
            for i in range(20):
                for j in range(i):
                    value=rdMolAlign.GetBestAlignmentTransform(molecules[i],molecules[j],map=maps,reflect=False,numThreads=1)[0]
                    matrix[i,j]=matrix[j,i]=value
            selected={str(c):separated_subset(matrix,c) for c in (.1,.25,.5)}
            rows.append(dict(case=tasks[0]['model_key'],campaign=str(path.parent),declared_starts=20,
                complete_source_search=plan.get('source_search_complete',True),heavy_atoms=template.GetNumAtoms(),
                stereochemistry_preserving_automorphisms=len(matches),heavy_atom_RMSD_A=matrix.tolist(),
                energy_ordered_separated_input_indices=selected,
                separated_input_counts={k:len(v) for k,v in selected.items()},
                DFT_minima_count_established=False,thermodynamic_convergence_established=False))
            sources[str(path)]=api.digest(path)
            atomic_write_json(args.out,dict(status='input_geometry_diversity_audit_in_progress',cases=rows,source_sha256=sources))
    api.verify_sources(sources)
    report=dict(status='declared_conformer_input_geometry_diversity_checked',checked=time.time(),cases=rows,
        source_sha256=sources,rdkit_version=rdBase.rdkitVersion,quantum_jobs_started=0,
        scope='Greedy energy-ordered separated input subsets, not an exact cluster count, DFT minima, probability or convergence certificate')
    atomic_write_json(args.out,report)
    print(json.dumps(dict(status=report['status'],cases=len(rows),counts_at_half_A={r['case']:r['separated_input_counts']['0.5'] for r in rows})))


if __name__=='__main__':main()
