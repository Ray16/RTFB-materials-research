"""Geometry-only descriptions and starting ensembles for pendant pi-core contacts.

These metrics classify geometry; they do not establish an interaction free energy,
an aggregate population, or loss of electrochemical accessibility.
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from rdkit import Chem, DistanceGeometry
from rdkit.Chem import AllChem, rdDistGeom


def enantiomer_representatives(smiles):
    """One representative per enantiomer pair, for achiral-pendant model controls."""
    found = {}
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        mirror = Chem.Mol(mol)
        for atom in mirror.GetAtoms():
            atom.InvertChirality()
        inverse = Chem.MolToSmiles(mirror, isomericSmiles=True)
        key = min(smi, inverse)
        found.setdefault(key, smi)
    return sorted(found.values())


def pi_core_groups(mol, fragment):
    """Map conjugated motif rings, excluding the polymer's phenyl rings.

    Ring membership is derived from the motif graph. This includes nonaromatic
    sp2 quinone rings, which an aromatic-atom-only filter would miss.
    """
    template = Chem.MolFromSmiles(fragment)
    query = Chem.MolFromSmarts(fragment)
    if template is None or query is None:
        raise ValueError('Invalid motif')
    ring_atoms = set()
    for ring in template.GetRingInfo().AtomRings():
        if len(ring) >= 5 and all(template.GetAtomWithIdx(i).GetHybridization() == Chem.HybridizationType.SP2 for i in ring):
            ring_atoms.update(ring)
    if len(ring_atoms) < 5:
        raise ValueError('No reviewed conjugated ring core')
    groups = {}
    for match in mol.GetSubstructMatches(query, uniquify=False, maxMatches=10000):
        ordered = tuple(match[i] for i in sorted(ring_atoms))
        groups.setdefault(frozenset(ordered), ordered)
    values = sorted(groups.values(), key=lambda g: min(g))
    if any(set(a)&set(b) for a,b in itertools.combinations(values,2)):
        raise ValueError('Ambiguous overlapping pi-core assignments')
    return [list(g) for g in values]


def plane(coordinates):
    xyz = np.asarray(coordinates, dtype=float)
    center = xyz.mean(axis=0)
    _, _, axes = np.linalg.svd(xyz-center, full_matrices=False)
    normal = axes[-1]
    rms = float(np.sqrt(np.mean(((xyz-center)@normal)**2)))
    return center, normal, axes[0], rms


def contact_metrics(positions, groups):
    xyz = np.asarray(positions, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or not np.isfinite(xyz).all():
        raise ValueError('Finite N by 3 coordinates required')
    output = []
    for (i,a),(j,b) in itertools.combinations(enumerate(groups),2):
        ca,na,_,ra = plane(xyz[a]); cb,nb,_,rb = plane(xyz[b])
        delta = cb-ca
        distance = float(np.linalg.norm(delta))
        angle = math.degrees(math.acos(float(np.clip(abs(na@nb),0,1))))
        separation = float((abs(delta@na)+abs(delta@nb))/2)
        slip = float(np.sqrt(max(0,distance**2-separation**2)))
        nearest = float(np.linalg.norm(xyz[a][:,None,:]-xyz[b][None,:,:],axis=2).min())
        # A transparent initial geometric heuristic, not a universal physical cutoff.
        stacked = angle <= 30 and 2.8 <= separation <= 4.2 and slip <= 3.5 and max(ra,rb) <= .2
        output.append(dict(site_i=i+1,site_j=j+1,centroid_distance_A=distance,
                           interplane_angle_deg=angle,plane_separation_A=separation,
                           lateral_offset_A=slip,nearest_core_atom_A=nearest,
                           core_planarity_rms_A=[ra,rb],stack_geometry=bool(stacked)))
    return output


def _constraint_bounds(mol, groups, mode, reference_xyz):
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(mol)
    if mode == 'unbiased' or len(groups) < 2:
        return bounds, []
    constraints = []
    if mode == 'open':
        for a,b in itertools.combinations(groups,2):
            for x,y in itertools.product(a,b):
                lo,hi=sorted((x,y))
                bounds[hi,lo]=max(bounds[hi,lo],5.5)
        if not DistanceGeometry.DoTriangleSmoothing(bounds):
            raise ValueError('Open bounds inconsistent with this backbone')
        return bounds,constraints
    # Three noncollinear core atoms constrain relative planes without distorting
    # the core itself. These biases are removed before xTB/DFT relaxation.
    core=np.asarray(reference_xyz)[groups[0]]
    center,normal,axis,_=plane(core)
    points=core-center
    pair=max(itertools.combinations(range(len(core)),2),key=lambda ij:np.linalg.norm(points[ij[0]]-points[ij[1]]))
    line=points[pair[1]]-points[pair[0]]
    third=max((i for i in range(len(core)) if i not in pair),key=lambda i:np.linalg.norm(np.cross(points[i]-points[pair[0]],line)))
    chosen=[pair[0],pair[1],third]
    ideal=[]
    for site in range(len(groups)):
        rotated=points if mode!='antiparallel' or site%2==0 else -points
        slip=site*2.0 if mode=='slipped' else 0.
        ideal.append(rotated+site*3.6*normal+slip*axis)
    for i,j in itertools.combinations(range(len(groups)),2):
        for ai,bi in itertools.product(chosen,chosen):
            a,b=groups[i][ai],groups[j][bi]
            target=float(np.linalg.norm(ideal[i][ai]-ideal[j][bi]))
            lower,upper=max(2.5,target-.45),target+.45
            lo,hi=sorted((a,b))
            bounds[hi,lo]=lower;bounds[lo,hi]=upper
            constraints.append((a,b,lower,upper))
    if not DistanceGeometry.DoTriangleSmoothing(bounds):
        raise ValueError('Stack bounds inconsistent with this backbone')
    return bounds,constraints


def starting_ensemble(smiles, fragment, mode, count=6, seed=20261019):
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
    groups=pi_core_groups(mol,fragment)
    params=AllChem.ETKDGv3();params.randomSeed=seed;params.numThreads=1
    params.maxIterations=1000;params.pruneRmsThresh=.3
    initial=Chem.Mol(mol)
    if AllChem.EmbedMolecule(initial,params)<0:
        params.useRandomCoords=True
        if AllChem.EmbedMolecule(initial,params)<0:
            raise RuntimeError('Reference embedding failed, including random-coordinate fallback')
    bounds,constraints=_constraint_bounds(mol,groups,mode,initial.GetConformer().GetPositions())
    params.SetBoundsMat(bounds)
    conformers=list(AllChem.EmbedMultipleConfs(mol,numConfs=count,params=params))
    if not conformers:
        params.useRandomCoords=True
        conformers=list(AllChem.EmbedMultipleConfs(mol,numConfs=count,params=params))
        if not conformers:raise RuntimeError('No embedded conformer for '+mode)
    props=AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:raise ValueError('MMFF parameterization unavailable for control')
    ranked=[]
    for cid in conformers:
        force=AllChem.MMFFGetMoleculeForceField(mol,props,confId=cid)
        for a,b,lower,upper in constraints:
            force.MMFFAddDistanceConstraint(a,b,False,lower,upper,20.)
        if mode=='open':
            for a,b in itertools.combinations(groups,2):
                force.MMFFAddDistanceConstraint(a[0],b[0],False,8.,25.,20.)
        force.Initialize();force.Minimize(maxIts=500)
        ranked.append((float(force.CalcEnergy()),cid))
    ranked.sort()
    return mol,groups,ranked
