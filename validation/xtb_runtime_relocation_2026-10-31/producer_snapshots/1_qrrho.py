"""CPU thermodynamic kernel with explicit internal vibrational modes.

The validated xTB 6.7.1 CLI filters all frequencies below 1 cm^-1, including
genuine internal modes. Its installed thermodyn kernel accepts an explicit mode
list. This module uses that kernel after excluding rigid motion geometrically.
It does not establish DFT stationarity, derivative accuracy or spin ordering.
"""
from __future__ import annotations

import atexit
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import selectors
import subprocess
import sys
import tempfile
import threading

import numpy as np
from redox.qm.xtb_runtime import configured_library

LIBRARY = configured_library()
LIBRARY_SHA256 = 'c405797f34ed96470f23158406e6d8a8e102681d3c3cbbaea16692f514f22853'
AUTORCM = 219474.63067
AMUTOAU = 1.660539040e-27 / 9.10938356e-31
HARTREE_EV = 27.211386245988

PROTOCOL = dict(schema=1, model='projected internal-mode xTB6.7.1 qRRHO',
                projection='mass-weighted SVD internal subspace; symmetric full-precision Cartesian Hessian',
                masses='installed xTB NIST isotopic-average atomic masses',
                rigid_mode_exclusion='geometric subspace; no frequency threshold',
                nuclear_kernel='installed GNU Fortran thermodyn; explicit internal modes',
                kernel_isolation='dedicated CPU subprocess; no PySCF/D4 Python imports',
                library_sha256=LIBRARY_SHA256, frequency_scale=1., temperature_K=298.15,
                rotor_cutoffs_cm=[25, 50, 100], primary_rotor_cutoff_cm=50,
                imaginary_mode_policy='absolute value only within the existing -5 cm^-1 minimum gate',
                zero_internal_frequency_policy='refuse; no artificial floor',
                symmetry='xTB symmetry detection for the exact input geometry',
                pressure='same translation convention as installed xTB thermodyn')


@lru_cache(maxsize=1)
def _checked_library():
    if hashlib.sha256(LIBRARY.read_bytes()).hexdigest() != LIBRARY_SHA256:
        raise RuntimeError('xTB thermodynamic library differs from the qualified build')
    return str(LIBRARY)


_WORKER=None
_LOCK=threading.Lock()
_LAST_WORKER=None


def _close_worker():
    global _WORKER
    if _WORKER is None or _WORKER['owner']!=os.getpid():return
    process=_WORKER['process']
    try:
        try:process.stdin.close()
        except OSError:pass
        try:process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            process.terminate();process.wait(timeout=3)
    finally:
        process.stdout.close();_WORKER['stderr'].close();_WORKER=None


atexit.register(_close_worker)


def _request(action,**arguments):
    """Keep incompatible Fortran libraries out of the caller's loader namespace."""
    global _WORKER,_LAST_WORKER
    _checked_library()
    source=Path(__file__).with_name('xtb_kernel.py')
    expected_source=hashlib.sha256(source.read_bytes()).hexdigest()
    with _LOCK:
        if _WORKER is not None and _WORKER['owner']!=os.getpid():_WORKER=None
        if _WORKER is None:
            env=dict(os.environ)
            for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
                env[key]='1'
            env['LD_LIBRARY_PATH']=str(LIBRARY.parent)+os.pathsep+env.get('LD_LIBRARY_PATH','')
            env.pop('LD_PRELOAD',None)
            error=tempfile.TemporaryFile(mode='w+t')
            process=subprocess.Popen([sys.executable,'-S',str(source),'--library',str(LIBRARY),
                                      '--sha256',LIBRARY_SHA256],stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                                     stderr=error,text=True,bufsize=1,env=env)
            _WORKER=dict(owner=os.getpid(),process=process,stderr=error,sequence=0)
        worker=_WORKER;process=worker['process'];worker['sequence']+=1
        request=dict(id=worker['sequence'],action=action,**arguments)
        try:
            process.stdin.write(json.dumps(request,allow_nan=False)+'\n');process.stdin.flush()
            with selectors.DefaultSelector() as selector:
                selector.register(process.stdout,selectors.EVENT_READ)
                if not selector.select(timeout=30):raise RuntimeError('Thermodynamic CPU worker response timed out')
            line=process.stdout.readline()
            if not line:raise RuntimeError('Thermodynamic CPU worker stopped before responding')
            response=json.loads(line)
            if (response.get('id')!=request['id'] or response.get('status')!='ok'
                    or response['worker']['pid']!=process.pid
                    or response['worker']['source_sha256']!=expected_source
                    or response['worker']['library_sha256']!=LIBRARY_SHA256):
                raise RuntimeError('Thermodynamic CPU worker rejected the request: '+str(response.get('error')))
            _LAST_WORKER=response['worker']
            return response['result']
        except BaseException as exc:
            worker['stderr'].seek(0);details=worker['stderr'].read()[-2000:]
            _close_worker()
            raise RuntimeError(str(exc)+((': '+details) if details else '')) from exc


def kernel_execution():
    """Ownership of the CPU helper that evaluated the most recent partition."""
    return dict(_LAST_WORKER or {})


def atomic_masses(numbers):
    """Return the exact masses used by the qualified xTB build, in amu."""
    values = np.asarray(numbers)
    if values.ndim != 1 or any(int(z) != z or not 1 <= z <= 118 for z in values):
        raise ValueError('atomic numbers between 1 and 118 required')
    return np.array(_request('masses',numbers=[int(z) for z in values]))


def internal_vibrations(matrix, coordinates_A, masses_amu):
    """Diagonalize only internal motion; retain signed unstable frequencies."""
    xyz = np.asarray(coordinates_A, dtype=float)
    mass = np.asarray(masses_amu, dtype=float)
    hessian = np.asarray(matrix, dtype=float)
    natm = len(xyz)
    if (xyz.shape != (natm, 3) or mass.shape != (natm,) or natm < 1
            or hessian.shape != (3 * natm, 3 * natm) or np.any(mass <= 0)
            or not all(np.isfinite(v).all() for v in (xyz, mass, hessian))):
        raise ValueError('finite Cartesian Hessian/coordinates and positive masses required')
    if abs(hessian - hessian.T).max() > 1e-6:
        raise ValueError('input Hessian fails the existing symmetry gate')
    centered = xyz - np.average(xyz, axis=0, weights=mass)
    roots = np.sqrt(mass)
    vectors = []
    for axis in np.eye(3):
        vectors.append((roots[:, None] * axis).ravel())
        vectors.append((roots[:, None] * np.cross(axis, centered)).ravel())
    full, singular, _ = np.linalg.svd(np.column_stack(vectors), full_matrices=True)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-10))
    internal = full[:, rank:]
    repeated = np.repeat(np.sqrt(mass * AMUTOAU), 3)
    weighted = (hessian + hessian.T) / (2 * np.outer(repeated, repeated))
    eigenvalues = np.linalg.eigvalsh(internal.T @ weighted @ internal)
    frequencies = np.sign(eigenvalues) * np.sqrt(abs(eigenvalues)) * AUTORCM
    return frequencies, rank


def inertia_properties(coordinates_A, masses_amu):
    """The installed axis2 inertia convention, with input coordinates in Å."""
    coords, mass = np.asarray(coordinates_A), np.asarray(masses_amu)
    centered = coords - np.average(coords, axis=0, weights=mass)
    second = np.einsum('i,ij,ik->jk', mass, centered, centered)
    inertia = np.eye(3) * np.trace(second) - second
    # Replay axis2's tiny perturbation rather than silently changing its model.
    inertia += np.array([[1., 2., 4.], [2., 3., 5.], [4., 5., 6.]]) * 1e-10
    moments = np.linalg.eigvalsh(inertia)
    rotational = np.array([16.8576522 / v if v >= 3e-4 else 0. for v in moments])[::-1]
    average = float(1e-47 * 1.66053 * sum(moments) / 3)
    return rotational, average, bool(rotational[-1] < 1e-10)


def kernel_thermal(frequencies_cm, rotational_constants_cm, average_moment_SI,
                   total_mass_amu, symmetry_number, temperature_K=298.15,
                   cutoff_cm=50, *, linear=False, atom=False):
    """Evaluate the installed kernel without the CLI's 1 cm^-1 mode filter."""
    frequencies = np.asarray(frequencies_cm, dtype=float)
    if (frequencies.ndim != 1 or not np.isfinite(frequencies).all()
            or np.any(frequencies <= 0)):
        raise ValueError('explicit strictly positive internal frequencies required')
    parameters = [*rotational_constants_cm, average_moment_SI, total_mass_amu,
                  symmetry_number, temperature_K, cutoff_cm]
    if (len(parameters) != 8 or not np.isfinite(parameters).all()
            or min(total_mass_amu, symmetry_number, temperature_K, cutoff_cm) <= 0
            or int(symmetry_number) != symmetry_number):
        raise ValueError('finite inertia, positive integer symmetry and thermodynamic parameters required')
    result=_request('thermal',parameters=dict(frequencies_cm=frequencies.tolist(),
        rotational_constants_cm=[float(v) for v in rotational_constants_cm],
        average_moment_SI=float(average_moment_SI),total_mass_amu=float(total_mass_amu),
        symmetry_number=int(symmetry_number),temperature_K=float(temperature_K),cutoff_cm=float(cutoff_cm),
        linear=bool(linear),atom=bool(atom)))
    if not all(np.isfinite(v) for v in result.values()):
        raise RuntimeError('nonfinite thermodynamic-kernel output')
    return result


def nuclear_thermals(matrix, coordinates_A, numbers, symmetry_number):
    """Compute every internal mode at 25/50/100 cm^-1; apply no new physical gate."""
    masses = atomic_masses(numbers)
    frequencies, rank = internal_vibrations(matrix, coordinates_A, masses)
    if len(frequencies) and (min(frequencies) < -5 or np.any(frequencies == 0)):
        raise ValueError('unstable or exactly zero internal mode; refine before thermochemistry')
    rotational, average, linear = inertia_properties(coordinates_A, masses)
    if len(numbers) > 1 and ((rank == 5) != linear):
        raise RuntimeError('geometric and inertia linearity classifications disagree')
    values = {str(cutoff): kernel_thermal(abs(frequencies), rotational, average, float(sum(masses)),
                                        symmetry_number, cutoff_cm=cutoff, linear=linear,
                                        atom=len(numbers) == 1)['G_Ha'] * HARTREE_EV
              for cutoff in PROTOCOL['rotor_cutoffs_cm']}
    return dict(frequencies_cm=frequencies.tolist(), rigid_rank=rank,
                internal_mode_count=len(frequencies), retained_soft_mode_indices=np.where(abs(frequencies) <= 1)[0].tolist(),
                imaginary_noise_mode_indices=np.where(frequencies < 0)[0].tolist(),
                nuclear_by_cutoff_eV=values, rotational_constants_cm=rotational.tolist(),
                average_moment_SI=average, atomic_masses_amu=masses.tolist())
