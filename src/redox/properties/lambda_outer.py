"""Outer-sphere (solvent) reorganization energy lambda_out — two continuum methods, reusable
for any molecule and any solvent (MeCN defaults from config/electrolyte.py).

Both are LINEAR-RESPONSE CONTINUUM methods and give the 1-body (electrochemical / half-reaction)
value. For homogeneous self-exchange use n_body=2 (== 2x at infinite separation; a contact cavity
would reduce this — see outer_sphere_validation.md). They are VERTICAL and use the fast/slow
(Pekar) partition, so they are strictly orthogonal to our GAS-phase 4-point lambda_in (no
double-counting). Optical dielectric = n**2 is used explicitly.

Methods
-------
1. born_lambda_outer : single-sphere Born/Marcus,
       lambda_o = z^2 e^2 / (8 pi eps0 a) * (1/eps_op - 1/eps_s)
   `a` = effective (SASA) radius. Cheap, shape-blind. (Same as solvated_reorg.lambda_outer_eV.)

2. pcm_lambda_outer  : MOLECULAR-CAVITY nonequilibrium continuum. Uses the REAL cavity + REAL
   charge-transfer potential. From the PCM apparent-charge operators q(eps)=K(eps)^-1 R(eps) v,
   the slow-polarization (Pekar) reorganization for the charge-transfer surface potential
   dV = V(ox) - V(red) on the cavity is
       lambda_o = 1/2 * dV . [ q_s(dV) - q_op(dV) ]
   with q_s built at eps_s and q_op at eps_op. Reduces to Born for a spherical cavity (the
   limiting-case test `sphere_sanity_check`: unit point charge in a single-atom cavity,
   reproduces analytic Born to <1e-5 relative). This is the method that tests whether the
   Born SPHERE is adequate for non-spherical ions.

Neither closes the continuum-vs-explicit gap (see outer_sphere_validation.md); they are the two
CONTINUUM points. Validate against explicit/experiment separately.
"""
from __future__ import annotations
import math

_E = 1.602176634e-19        # C
_EPS0 = 8.8541878128e-12    # F/m
_BOHR = 0.52917721067       # Angstrom
HARTREE_EV = 27.211386245988


# ----------------------------------------------------------------------------- solvent config
def solvent_constants(name: str | None = None) -> dict:
    """(eps_s, eps_op, refractive_index, name) for the electrolyte, from config/electrolyte.py
    (config/project.json). No silent fallback: a missing/broken config raises."""
    from redox.core.common import load_config
    S = load_config("electrolyte").SOLVENT
    return dict(eps_s=float(S.get("eps_r", S.get("eps"))), eps_op=float(S["eps_optical"]),
                n=float(S.get("refractive_index", S["eps_optical"] ** 0.5)),
                name=S.get("name", "config"))


# common solvents for validation across media (eps_s, refractive index n)
SOLVENTS = {
    "acetonitrile": (37.5, 1.344),
    "water":        (78.36, 1.333),
    "thf":          (7.58, 1.407),
    "dcm":          (8.93, 1.424),
    "dmf":          (36.7, 1.430),
    "dmso":         (46.7, 1.479),
}


# ----------------------------------------------------------------------------- 1. Born sphere
def born_lambda_outer(radius_A: float, eps_s: float, eps_op: float, z: int = 1,
                      n_body: int = 1) -> float:
    """Single-sphere Born/Marcus outer-sphere lambda (eV). n_body=2 -> self-exchange (x2)."""
    a = radius_A * 1e-10
    pekar = (1.0 / eps_op) - (1.0 / eps_s)
    lam = (z * z) * _E * pekar / (8.0 * math.pi * _EPS0 * a)   # eV (single-ion / 1-body)
    return n_body * lam


def born_two_sphere(a1_A: float, a2_A: float, d_A: float, eps_s: float, eps_op: float,
                    z: int = 1) -> float:
    """Marcus TWO-SPHERE outer-sphere lambda (eV) for a bimolecular (self-exchange) ET:
        lambda_o = z^2 e^2/(4 pi eps0) * (1/(2 a1) + 1/(2 a2) - 1/d) * (1/eps_op - 1/eps_s)
    a1,a2 = donor/acceptor radii (A); d = center-to-center distance (A).
    Limits (self-exchange a1=a2=a):  d=2a (contact) -> equals single-sphere Born (1-body);
    d->inf -> 2x single-sphere. Requires d >= a1+a2 (no cavity overlap)."""
    a1, a2, d = a1_A * 1e-10, a2_A * 1e-10, d_A * 1e-10
    pekar = (1.0 / eps_op) - (1.0 / eps_s)
    geom = 1.0 / (2.0 * a1) + 1.0 / (2.0 * a2) - 1.0 / d
    return (z * z) * _E * pekar / (4.0 * math.pi * _EPS0) * geom


def sasa_radius_A(smiles: str, probe_radius_A: float = 1.3, seed: int = 0xC0FFEE) -> float | None:
    """Solvent-accessible-surface equivalent-sphere radius (Angstrom)."""
    from redox.properties.solvated_reorg import hydrodynamic_radius_A
    return hydrodynamic_radius_A(smiles, probe_radius_A=probe_radius_A, method="sasa", seed=seed)


# ------------------------------------------------------ 2. molecular-cavity nonequilibrium PCM
def _pcm_operators(mol, eps: float):
    """Build a PySCF PCM at dielectric `eps`; return (cm, K, R, v_grids_n)."""
    from pyscf.solvent import pcm as _pcm
    cm = _pcm.PCM(mol)
    cm.eps = float(eps)
    cm.method = "IEF-PCM"
    cm.build()
    K = cm._intermediates["K"]; R = cm._intermediates["R"]
    return cm, K, R, cm.v_grids_n


def _surface_potential(cm, dm):
    """Solute electrostatic potential on the cavity grid: V = V_nuc - V_elec."""
    import numpy as np
    v_e = cm._get_v(np.asarray(dm).reshape(1, dm.shape[-1], dm.shape[-1]))[0]
    return cm.v_grids_n - v_e


def _q_sym(K, R, v):
    """PySCF apparent surface charges for surface potential v (symmetrized as in
    pyscf.solvent.pcm._get_vind). Convention: R = -f(eps)(...), so q is OPPOSITE in sign to v
    and the solvation energy 1/2 q.v is NEGATIVE (stabilizing)."""
    import numpy as np
    q = np.linalg.solve(K, R.dot(v))
    qt = R.T.dot(np.linalg.solve(K.T, v))
    return 0.5 * (q + qt)


def _pekar_lambda_ha(K_s, R_s, K_op, R_op, dV):
    """Nonequilibrium (Pekar) reorganization for a charge-transfer surface potential dV:
        lambda_o = E_solv(eps_op) - E_solv(eps_s) = 1/2 dV.(q_op - q_s)   (Hartree, >= 0)
    i.e. the slow-polarization energy that must be paid. NO abs(): a negative value means a
    sign/operator error and is raised, not hidden."""
    lam = 0.5 * float(dV.dot(_q_sym(K_op, R_op, dV) - _q_sym(K_s, R_s, dV)))
    if lam < -1e-8:
        raise ValueError(f"negative Pekar lambda_o ({lam:.3e} Ha): PCM sign/operator error")
    return lam


def pcm_lambda_outer(mol, dm_red, dm_ox, eps_s: float, eps_op: float, n_body: int = 1) -> float:
    """Molecular-cavity nonequilibrium (Pekar) outer-sphere lambda (eV), 1-body by default.

    dm_red, dm_ox : converged density matrices of the two redox states at the SAME geometry
                    (total DM; for UKS pass dm_a+dm_b). Uses the charge-transfer surface
                    potential dV = V(ox) - V(red) and the slow (eps_s minus eps_op) response.
    """
    import numpy as np

    def total(dm):
        dm = np.asarray(dm)
        return dm[0] + dm[1] if dm.ndim == 3 else dm

    dm_r, dm_o = total(dm_red), total(dm_ox)
    cm_s, K_s, R_s, _ = _pcm_operators(mol, eps_s)          # cavity is geometry-only
    dV = _surface_potential(cm_s, dm_o) - _surface_potential(cm_s, dm_r)
    _, K_op, R_op, _ = _pcm_operators(mol, eps_op)
    return n_body * _pekar_lambda_ha(K_s, R_s, K_op, R_op, dV) * HARTREE_EV


def sphere_sanity_check(eps_s: float = 37.5, eps_op: float = 1.806, element: str = "Ne"):
    """REAL limiting-case test of the PCM operator path: a unit POINT charge at the centre of a
    single-atom (spherical) cavity. The Pekar lambda from the same K/R operators and the same
    _pekar_lambda_ha used for molecules must reproduce analytic Born for that cavity radius.
    Returns dict(born_eV, pcm_eV, radius_A, rel_err)."""
    import numpy as np
    from pyscf import gto
    mol = gto.M(atom=f"{element} 0 0 0", basis="sto-3g", verbose=0)
    cm_s, K_s, R_s, vn = _pcm_operators(mol, eps_s)
    _, K_op, R_op, _ = _pcm_operators(mol, eps_op)
    dV = vn / float(mol.atom_charges()[0])                   # potential of a +1 point charge
    lam = _pekar_lambda_ha(K_s, R_s, K_op, R_op, dV) * HARTREE_EV
    r_A = float(np.linalg.norm(cm_s.surface["grid_coords"], axis=1).mean()) * _BOHR
    born = born_lambda_outer(r_A, eps_s, eps_op)
    return dict(born_eV=round(born, 5), pcm_eV=round(lam, 5), radius_A=round(r_A, 4),
                rel_err=round((lam - born) / born, 5))
