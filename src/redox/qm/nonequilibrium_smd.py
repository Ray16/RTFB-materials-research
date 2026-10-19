"""Experimental frozen-slow extension of native SMD electrostatics.

The SMD apparent-charge response is IEF-PCM on SMD element/solvent radii.
Its separate CDS term is independent of the electronic density at fixed
coordinates. This extension retains the native SMD cavity and CDS functional;
only electrostatic inertial polarization is frozen. CDS is not an additional
outer-sphere lambda. Its treatment away from equilibrium is an explicit model
assumption, not a validated model of solvent structural dynamics.

Native solvent descriptors determine eps_static and eps_optical=n**2 unless
explicitly overridden; these must be recorded because project-level constants
can differ from the SMD library. No production promotion is done here.
"""
from __future__ import annotations
import hashlib
import types

import numpy as np

from redox.qm.nonequilibrium_pcm import host_array, surface_potential


def smd_object(mol, solvent="acetonitrile", backend="cpu", eps=None, order=41):
    if backend=="cpu":
        from pyscf.solvent.smd import SMD
    elif backend=="gpu":
        from gpu4pyscf.solvent.smd import SMD
    else:
        raise ValueError("backend must be cpu or gpu")
    cm=SMD(mol,solvent=solvent)
    cm.lebedev_order=order
    if eps is not None:
        if not np.isfinite(eps) or eps<1:
            raise ValueError("finite dielectric >=1 required")
        cm.eps=float(eps)
    cm.build()
    return cm


def smd_dielectrics(solvent="acetonitrile"):
    from pyscf.solvent.smd import solvent_db
    desc=solvent_db[solvent]
    return dict(eps_s=float(desc[5]),eps_op=float(desc[0]**2),
                refractive_index=float(desc[0]),source="native_pyscf_smd_solvent_db")


def attach_frozen_slow_smd(mf, reference_dm, solvent="acetonitrile",backend="cpu",
                            eps_s=None,eps_op=None,order=41):
    if reference_dm.shape[-2:]!=(mf.mol.nao,mf.mol.nao):
        raise ValueError("reference density is not in the target AO basis")
    native_constants=smd_dielectrics(solvent)
    eps_s=native_constants['eps_s'] if eps_s is None else float(eps_s)
    eps_op=native_constants['eps_op'] if eps_op is None else float(eps_op)
    if not np.isfinite([eps_s,eps_op]).all() or not 1<=eps_op<=eps_s:
        raise ValueError("require1<=optical<=static dielectric")
    static=smd_object(mf.mol,solvent,backend,eps_s,order)
    optical=smd_object(mf.mol,solvent,backend,eps_op,order)
    for key in ('grid_coords','charge_exp'):
        if not np.array_equal(host_array(static.surface[key]),host_array(optical.surface[key])):
            raise RuntimeError("static and optical SMD cavities differ")
    static._get_vind(reference_dm)
    optical._get_vind(reference_dm)
    slow=(static._intermediates['q_sym']-optical._intermediates['q_sym']).copy()
    vslow=optical._get_vmat(slow.reshape(1,-1))[0]
    self_work=.5*float(host_array(surface_potential(optical,reference_dm)@slow))
    native=optical._get_vind

    def get_vind(self,dm):
        energy,matrix=native(dm)
        interaction=surface_potential(self,dm)@slow
        return energy+interaction-self_work,matrix+vslow

    optical._get_vind=types.MethodType(get_vind,optical)
    # Native .SMD adds CDS exactly once in energy_elec. _get_vind itself
    # returns electrostatic energy only, exactly as the equilibrium class.
    mf=mf.SMD(optical)
    mf.nonequilibrium_model=dict(partition='Pekar',model='SMD_cavity_frozen_electrostatics',
        solvent=solvent,eps_s=eps_s,eps_op=eps_op,sasa_ng=optical.sasa_ng,
        cds_assumption='native_density_independent_geometry_functional',
        reference_dm_sha256=hashlib.sha256(host_array(reference_dm).tobytes()).hexdigest(),
        validated_physical_accuracy=False)
    return mf


def cds_relaxation_difference(initial_smd,final_smd):
    """CDS contribution to final-state relaxation from initial to final coords.

    Equivalent directional totals have opposite CDS differences, so their
    symmetric mean cancels CDS exactly. It still affects optimized geometries.
    """
    if initial_smd.solvent!=final_smd.solvent:
        raise ValueError("cannot compare CDS from different solvents")
    return float(initial_smd.get_cds())-float(final_smd.get_cds())


def independent_hcore_smd_reference(mf,reference_dm,solvent='acetonitrile',order=41):
    """CPU reference assembly preserving native SMD+CDS without kernel override."""
    constants=smd_dielectrics(solvent)
    static=smd_object(mf.mol,solvent,'cpu',constants['eps_s'],order)
    optical=smd_object(mf.mol,solvent,'cpu',constants['eps_op'],order)
    static._get_vind(reference_dm);optical._get_vind(reference_dm)
    slow=static._intermediates['q_sym']-optical._intermediates['q_sym']
    vslow=optical._get_vmat(slow.reshape(1,-1))[0]
    self_work=.5*float(surface_potential(optical,reference_dm)@slow)
    native_hcore=mf.get_hcore
    def get_hcore(self,mol=None):return native_hcore(mol)+vslow
    mf.get_hcore=types.MethodType(get_hcore,mf)
    return mf.SMD(optical),float(optical.v_grids_n@slow)-self_work
