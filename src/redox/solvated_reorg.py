"""Solvated (total) reorganization energy = inner-sphere + outer-sphere.

    lambda_total = lambda_i + lambda_o                              [Sharma 2021, Eq. 9]

- lambda_i (inner-sphere): the geometric/vibrational term we already compute by the Nelsen
  4-point scheme (src/redox/reorg.py), GAS-phase. This is the "gas phase" reorganization.
- lambda_o (outer-sphere): the solvent (Born/Marcus continuum) term the inner-sphere misses.
  From Sharma et al., PNAS 2021 (Electrochemical implications of modulating the solvation
  shell...), Eq. 12-13 -- the Marcus/Born expression with the Pekar factor (1/eps_op - 1/eps_r):

    lambda_o = (z^2 e^2) / (8 pi eps0 a) * (1/eps_op - 1/eps_r)     [Eq. 12/13]

  z  = electrons transferred (=1 for a 1e couple); e = elementary charge; eps0 = vacuum
  permittivity; a = effective radius of the species (Sharma uses the hydrodynamic radius; we
  use the effective radius from the molecular volume, a = (3V/4pi)^(1/3), as a computable
  proxy -- absolute lambda_o is therefore approximate, but the trends between species and the
  gas->solvated shift are robust). eps_op = optical permittivity (= n^2), eps_r = static
  permittivity. Sharma's numbers are for WATER; we read MeCN from config/electrolyte.py
  (eps_r=37.5, eps_op=1.806).

The "solvated reorganization energy" for the 2x2 SA-vs-lambda landscape is lambda_total; the
"gas phase" value is lambda_i alone.
"""
from __future__ import annotations
import math

_E = 1.602176634e-19          # elementary charge, C
_EPS0 = 8.8541878128e-12      # vacuum permittivity, F/m


def lambda_outer_eV(radius_A: float, eps_op: float, eps_r: float, z: int = 1) -> float:
    """Outer-sphere (solvent) reorganization energy in eV, Sharma Eq. 12/13.

    radius_A : effective species radius in Angstrom (see effective_radius_A).
    eps_op   : optical permittivity (n^2). eps_r: static permittivity. z: electrons (1e -> 1).
    Returns lambda_o in eV (>= 0; the Pekar factor 1/eps_op - 1/eps_r > 0 since eps_op < eps_r).
    """
    a = radius_A * 1e-10                       # Angstrom -> m
    pekar = (1.0 / eps_op) - (1.0 / eps_r)
    # e^2/(8 pi eps0 a) * pekar is in Joules; dividing by e (once) gives eV, leaving one e:
    return (z * z) * _E * pekar / (8.0 * math.pi * _EPS0 * a)


def _embed(smiles: str, seed: int = 0xC0FFEE):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    if AllChem.EmbedMolecule(m, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(m, p) != 0:
            return None
    return m


def hydrodynamic_radius_A(smiles: str, probe_radius_A: float, method: str = "sasa",
                          seed: int = 0xC0FFEE) -> float | None:
    """Estimate the hydrodynamic radius r_hyd (Angstrom) for Sharma Eq. 13.

    method="sasa" (DEFAULT, preferred): roll a solvent-sized probe over the solute and take
    the equivalent-sphere radius of its SOLVENT-ACCESSIBLE surface, r = sqrt(SASA / 4pi). This
    includes the solvent contact shell, so it is a far better r_hyd proxy than the bare
    molecular volume (which omits the shell and underestimates r_hyd -> overestimates lambda_o).
    method="volume": bare van-der-Waals equivalent-sphere radius a = (3V/4pi)^(1/3) (kept for
    comparison / fallback; ignores probe_radius_A).

    Caveat: still an ESTIMATE. The true r_hyd is defined by diffusion (Stokes-Einstein
    r_hyd = kT/6*pi*eta*D); without a measured or MD diffusion coefficient, the SASA radius is
    the best purely-structural proxy. Returns None if embedding/SASA fails."""
    m = _embed(smiles, seed)
    if m is None:
        return None
    if method == "volume":
        from rdkit.Chem import AllChem
        vol = AllChem.ComputeMolVolume(m)
        return (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0)
    # SASA path
    try:
        from rdkit.Chem import rdFreeSASA
        radii = rdFreeSASA.classifyAtoms(m)
        opts = rdFreeSASA.SASAOpts(); opts.probeRadius = float(probe_radius_A)
        sasa = rdFreeSASA.CalcSASA(m, radii, opts=opts)          # Angstrom^2
        return math.sqrt(sasa / (4.0 * math.pi))
    except Exception:
        from rdkit.Chem import AllChem
        vol = AllChem.ComputeMolVolume(m)
        return (3.0 * vol / (4.0 * math.pi)) ** (1.0 / 3.0)


# Back-compat alias: previous callers used effective_radius_A (bare vdW volume radius).
def effective_radius_A(smiles: str, seed: int = 0xC0FFEE) -> float | None:
    """DEPRECATED name: bare vdW volume radius. Prefer hydrodynamic_radius_A(method='sasa')."""
    return hydrodynamic_radius_A(smiles, probe_radius_A=0.0, method="volume", seed=seed)


def lambda_total_eV(lambda_i_eV: float, smiles: str, eps_op: float, eps_r: float,
                    probe_radius_A: float, z: int = 1, radius_method: str = "sasa") -> dict:
    """Total solvated reorganization energy lambda_i + lambda_o for a species.
    Uses the SASA-based hydrodynamic radius (probe_radius_A) by default. Returns
    dict(r_hyd_A, radius_method, lambda_o_eV, lambda_total_eV); lambda_o=None if radius failed."""
    a = hydrodynamic_radius_A(smiles, probe_radius_A, method=radius_method)
    if a is None:
        return dict(r_hyd_A=None, radius_method=radius_method,
                    lambda_o_eV=None, lambda_total_eV=None)
    lo = lambda_outer_eV(a, eps_op, eps_r, z)
    return dict(r_hyd_A=round(a, 3), radius_method=radius_method,
                lambda_o_eV=round(lo, 4), lambda_total_eV=round(lambda_i_eV + lo, 4))
