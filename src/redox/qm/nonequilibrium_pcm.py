"""Diagnostic state-specific nonequilibrium PCM for vertical electron transfer.

Pekar partition, You et al., JCP 143, 204104 (2015), Eqs. 16, 21--23:
    G[D;D0] = E_vac[D] + .5 V Qop V + V (Qs-Qop)V0
              - .5 V0 (Qs-Qop)V0.
V includes nuclei and electrons. Q is symmetrized, exactly as in native
PySCF's variational PCM energy. Optical polarization relaxes in SCF; the
slow reaction field and its polarization work remain fixed. Freezing the
entire native PCM potential would incorrectly freeze the optical response.

This is electrostatic IEF-PCM, not equilibrium SMD and not a gas inner-only
energy. No gradients, geometry optimization, cavity transport or production
ranking promotion are implemented. Every vertical event MUST use the same
coordinates, basis and cavity for its initial and final electronic states.
"""
from __future__ import annotations

import numpy as np


def host_array(x):
    return np.asarray(x.get() if hasattr(x, "get") else x)


def total_density(dm):
    if dm.ndim == 3:
        if dm.shape[0] != 2:
            raise ValueError("expected RKS or alpha/beta UKS density")
        return dm[0] + dm[1]
    return dm


def surface_potential(pcm, dm):
    dm = total_density(dm)
    if hasattr(pcm, "_get_vgrids"):
        return pcm._get_vgrids(dm, with_nuc=True)[0]
    return pcm.v_grids_n - pcm._get_v(dm.reshape(1, *dm.shape))[0]


def pcm_object(mol, eps, backend="cpu", method="IEF-PCM", order=29):
    if not 1 <= float(eps) < float("inf"):
        raise ValueError("finite dielectric >=1 required")
    if backend == "gpu":
        from gpu4pyscf.solvent.pcm import PCM
    elif backend == "cpu":
        from pyscf.solvent.pcm import PCM
    else:
        raise ValueError("backend must be cpu or gpu")
    cm = PCM(mol)
    cm.eps = float(eps)
    cm.method = method
    cm.lebedev_order = order
    cm.build()
    return cm


def frozen_slow_terms(mol, reference_dm, eps_s, eps_op, backend="cpu",
                      method="IEF-PCM", order=29):
    """Return optical PCM, slow charges/matrix and invariant self-work.

    A density matrix alone has no geometry label; callers are responsible for
    passing it in the exact AO basis and coordinates of `mol`. Shape checks
    reject obvious mismatches and cavity checks reject differing grids.
    """
    if not 1 <= eps_op <= eps_s or not np.isfinite([eps_s, eps_op]).all():
        raise ValueError("require 1 <= optical <= static dielectric")
    if reference_dm.shape[-2:] != (mol.nao, mol.nao):
        raise ValueError("reference density is not in the target AO basis")
    static = pcm_object(mol, eps_s, backend, method, order)
    optical = pcm_object(mol, eps_op, backend, method, order)
    for key in ("grid_coords", "charge_exp"):
        if not np.array_equal(host_array(static.surface[key]),
                              host_array(optical.surface[key])):
            raise ValueError("static/optical PCM cavities differ")
    static._get_vind(reference_dm)
    optical._get_vind(reference_dm)
    slow_q = (static._intermediates["q_sym"]
              - optical._intermediates["q_sym"]).copy()
    slow_vmat = optical._get_vmat(slow_q.reshape(1, -1))[0]
    ref_v = surface_potential(optical, reference_dm)
    self_work = float(host_array(ref_v @ slow_q)) * .5
    return optical, slow_q, slow_vmat, self_work


def attach_frozen_slow(mf, reference_dm, eps_s, eps_op, backend="cpu",
                       method="IEF-PCM", order=29):
    """Attach optical PCM with fixed slow field and correct energy functional.

    Native PCM.kernel computes the optical energy and its AO derivative. The
    fixed slow field contributes its *full* interaction, not half, less the
    constant work needed to form the original slow polarization.
    """
    import types
    cm, qslow, vslow, self_work = frozen_slow_terms(
        mf.mol, reference_dm, eps_s, eps_op, backend, method, order)
    native = cm._get_vind

    def get_vind(self, dm):
        energy, matrix = native(dm)
        potential = surface_potential(self, dm)
        interaction = potential @ qslow
        return energy + interaction - self_work, matrix + vslow

    cm._get_vind = types.MethodType(get_vind, cm)
    cm.nonequilibrium_metadata = dict(
        partition="Pekar", eps_s=eps_s, eps_op=eps_op,
        model=method, lebedev_order=order, slow_self_work_Ha=self_work,
        reference_density_sha256=__import__("hashlib").sha256(
            host_array(reference_dm).tobytes()).hexdigest())
    # Passing dm to native .PCM would freeze the ENTIRE solvent; do not do so.
    return mf.PCM(cm)


def independent_hcore_reference(mf, reference_dm, eps_s, eps_op,
                               method="IEF-PCM", order=29):
    """CPU reference using native optical PCM plus a fixed external AO field.

    This path does not override PCM kernel. The slow electronic interaction
    is included in hcore, while the slow nuclear interaction and self-work
    are added separately. It provides an independent energy-bookkeeping/SCF
    check on attach_frozen_slow (same continuum physics, different assembly).
    """
    import types
    cm, qslow, vslow, self_work = frozen_slow_terms(
        mf.mol, reference_dm, eps_s, eps_op, "cpu", method, order)
    native_hcore = mf.get_hcore

    def get_hcore(self, mol=None):
        return native_hcore(mol) + vslow

    mf.get_hcore = types.MethodType(get_hcore, mf)
    mf = mf.PCM(cm)
    constant = float(cm.v_grids_n @ qslow) - self_work
    return mf, constant
