"""Explicitly experimental smooth-CDS context, isolated from primary records.

Electrostatic SMD, SCF and basis settings are unchanged. Only the non-electronic
surface-area term uses the installed experimental SWIG formulation. The model
identity differs from the legacy protocol, and no primary cache is written here.
"""
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

MODEL = dict(name='experimental_PySCF_SWIG_CDS_2030_v1', solvent='acetonitrile',
    electrostatic_model='SMD', electrostatic_surface_points=590,
    CDS_surface='PySCF experimental SWIG solvent-accessible area', CDS_points_per_sphere=2030,
    supported_elements=['C','H','N','O','Cl'], scientific_model_changed=True,
    primary_candidate_admission=False)


def component_evidence(root):
    path=Path(root)/'summary.json';row=json.loads(path.read_text())
    if row['status']!='experimental_smooth_CDS_component_study_complete':
        raise ValueError('Completed smooth-CDS component study required')
    api.verify_sources(row['source_sha256'])
    values=row['states']
    required={'contact','empty','capped','site2','site3','precursor','local','chloride','base','protonated','pair'}
    if not required.issubset(values):raise ValueError('Missing component qualification coverage')
    contact=values['contact']
    numeric=[contact['sampled_gradient']['maximum_error_Ha_Bohr'],
        contact['CPU_legacy_contact_reference_error_Ha'],contact['rigid_rotation_spread_kJ_mol'],
        *[r['finest_grid_change_kJ_mol'] for r in values.values()]]
    if not np.isfinite(numeric).all():
        raise ValueError('Missing or nonfinite component qualification evidence')
    if (contact['sampled_gradient']['tested_components']!=4
            or contact['sampled_gradient']['maximum_error_Ha_Bohr']>1e-6
            or contact['CPU_legacy_contact_reference_error_Ha']>1e-11
            or contact['rigid_rotation_spread_kJ_mol']>.5
            or max(r['finest_grid_change_kJ_mol'] for r in values.values())>.5):
        raise ValueError('Gradient, reference, quadrature or rotation component check failed')
    curvature=list(contact['contact_scan']['centered_curvature_Ha_Bohr2'].values())
    if not np.isfinite(curvature).all() or max(curvature)-min(curvature)>1e-5:
        raise ValueError('Smooth contact curvature has not stabilized')
    review=Path(root)/'independent_review.json'
    checked=json.loads(review.read_text())
    if checked['status']!='smooth_CDS_fixed_geometry_component_review_passed':
        raise ValueError('Independent library-bound component review required')
    api.verify_sources(checked['source_sha256'])
    return {str(path):api.digest(path),str(review):api.digest(review),
        **row['source_sha256'],**checked['source_sha256']}


class LazyCDS:
    def __init__(self, obj, energy, gradient, statistics):
        self.obj,self.energy,self.gradient,self.statistics=obj,energy,gradient,statistics
        self.values={}

    def __getitem__(self,index):
        if index not in (0,1):raise IndexError(index)
        if index not in self.values:
            value=self.energy(self.obj) if index==0 else self.gradient(self.obj)
            if not np.isfinite(value).all():raise ValueError('Nonfinite smooth-CDS value')
            self.values[index]=value
            self.statistics['energy_calls' if index==0 else 'gradient_calls']+=1
        return self.values[index]

    def __iter__(self):
        yield self[0]
        yield self[1]


@contextmanager
def smooth_cds(directory,component_root):
    from pyscf import gto
    from pyscf.solvent import smd_experiment as smooth
    from gpu4pyscf.solvent import smd

    evidence=component_evidence(component_root)
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=False)
    loader=api.load_module('private_experimental_CDS_gradient',ROOT/'scripts/validation/qualify_smooth_cds.py')
    derivative,gradient_source=loader.gradient_functions(smooth,directory/'experimental_gradient_source.py')
    profile=dict(MODEL,implementation_sha256={str(p):api.digest(p) for p in
        (Path(__file__).resolve(),Path(smooth.__file__),gradient_source)})
    profile_id=hashlib.sha256(json.dumps(profile,sort_keys=True).encode()).hexdigest()
    proof=dict(status='experimental_smooth_CDS_active',profile=profile,profile_id=profile_id,
        component_evidence=evidence,energy_calls=0,gradient_calls=0,primary_candidate_admission=False)
    original=smd.get_cds_legacy
    last_key,last_value=None,None
    def evaluate(solvent):
        nonlocal last_key,last_value
        if solvent.solvent.lower()!='acetonitrile':raise ValueError('Smooth pilot medium must be MeCN')
        if solvent.method.upper()!='SMD' or solvent.sasa_ng!=590:
            raise ValueError('Electrostatic SMD model/grid differs from this experimental profile')
        mol=solvent.mol
        descriptors=solvent.solvent_descriptors
        if descriptors is not None and not np.array_equal(descriptors,smooth.solvent_db['acetonitrile']):
            raise ValueError('Custom solvent descriptors differ from the qualified MeCN component model')
        symbols=tuple(mol.atom_symbol(i) for i in range(mol.natm))
        if set(symbols)-set(MODEL['supported_elements']):raise ValueError('Unsupported smooth-pilot element')
        coords=mol.atom_coords(unit='Bohr')
        coords=np.asarray(coords.get() if hasattr(coords,'get') else coords)
        key=(symbols,mol.charge,mol.spin,coords.tobytes())
        if key!=last_key:
            cpu=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',basis='sto-3g',
                charge=mol.charge,spin=mol.spin,verbose=0)
            obj=SimpleNamespace(mol=cpu,solvent='acetonitrile',_solvent='acetonitrile',
                solvent_descriptors=solvent.solvent_descriptors,sasa_ng=2030)
            last_key,last_value=key,LazyCDS(obj,smooth.get_cds,derivative,proof)
        return last_value
    smd.get_cds_legacy=evaluate
    try:
        yield proof
        proof['status']='experimental_smooth_CDS_context_finished_QC_separate'
    except BaseException as error:
        proof.update(status='experimental_smooth_CDS_failed_preserved',error=repr(error))
        raise
    finally:
        smd.get_cds_legacy=original
        atomic_write_json(directory/'model.json',proof)
