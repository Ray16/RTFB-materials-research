"""Source-addressed CPU postprocessing of already qualified native Hessians.

Original DFT records and their unprojected xTB logs remain immutable. A separate
record binds the corrected partition function to its geometry, raw native record,
Cartesian Hessian and qualified thermodynamic protocol. It never establishes
the backing DFT derivative or spin admission on its own.
"""
from __future__ import annotations

import fcntl
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

PRODUCER_SOURCE_BYTES=Path(__file__).read_bytes()

if __name__=='__main__':
    # The command owns its CPU allocation. Set pools before importing NumPy;
    # importing this module as a library does not change the caller's settings.
    preflight_parser=argparse.ArgumentParser(add_help=False)
    preflight_parser.add_argument('--threads',type=int,default=2)
    preflight_args,_=preflight_parser.parse_known_args()
    if preflight_args.threads<1:raise SystemExit('positive --threads required')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):
        os.environ[key]=str(preflight_args.threads)

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm import qrrho as Q
from redox.qm.vibrations import project_rigid_modes, write_xtb_hessian, parse_xtb_thermal_eV
QRRHO_SOURCE_BYTES=Path(Q.__file__).read_bytes()
KERNEL_SOURCE_BYTES=Path(Q.__file__).with_name('xtb_kernel.py').read_bytes()

CERTIFICATE = ROOT / 'config/projected_thermal_validation.json'
from redox.qm.xtb_runtime import configured_binary
XTB = configured_binary()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def certificate():
    proof = json.loads(CERTIFICATE.read_text())
    content_hash = proof.pop('certificate_content_sha256')
    actual = hashlib.sha256(json.dumps(proof, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if (actual != content_hash or proof.get('status') != 'independently_qualified'
            or proof.get('protocol') != Q.PROTOCOL
            or proof.get('qrrho_source_sha256') != sha(Path(Q.__file__))
            or proof.get('kernel_worker_source_sha256') != sha(Path(Q.__file__).with_name('xtb_kernel.py'))
            or proof.get('library_sha256') != sha(Q.LIBRARY)):
        raise RuntimeError('Projected thermodynamic implementation is not independently qualified')
    if proof.get('CLI_reproduction_bound_eV') != 1e-6:
        raise RuntimeError('Unqualified thermodynamic reproduction tolerance')
    for path, digest in proof['evidence_sha256'].items():
        if sha(ROOT / path) != digest:
            raise RuntimeError('Projected thermodynamic qualification evidence changed')
    return content_hash


def source_identity(native_path, geometry):
    native_path, geometry = Path(native_path), Path(geometry)
    raw = json.loads(native_path.read_text())
    hessian = native_path.with_suffix('') / 'hessian_cartesian.npy'
    ident = dict(geometry_sha256=sha(geometry), raw_native_sha256=sha(native_path),
                 Cartesian_Hessian_sha256=sha(hessian), charge=raw['identity']['charge'],
                 mult=raw['identity']['mult'], protocol=Q.PROTOCOL)
    if raw['identity']['geometry_sha256'] != ident['geometry_sha256']:
        raise RuntimeError('Native Hessian does not match the current geometry')
    return ident, raw, hessian


def cache_path(native_path, geometry):
    ident, _, _ = source_identity(native_path, geometry)
    key = hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16]
    return Path(native_path).parent.parent / 'projected_thermal' / (key + '.json'), ident


def _serialized_frequencies(path, masses):
    numbers = np.fromstring(' '.join(line for line in Path(path).read_text().splitlines()
                                     if not line.lstrip().startswith('$')), sep=' ')
    roots = np.repeat(np.sqrt(np.asarray(masses) * Q.AMUTOAU), 3)
    matrix = numbers.reshape(len(roots), len(roots)).T
    values = np.linalg.eigvalsh(matrix / np.outer(roots, roots))
    return np.sign(values) * np.sqrt(abs(values)) * Q.AUTORCM


def load_projected_thermal(native_path, geometry):
    """Read a corrected partition only; callers must also admit its raw Hessian."""
    try:
        proof = certificate()
        target, ident = cache_path(native_path, geometry)
        result = json.loads(target.read_text())
        if (result.get('status') != 'projected_internal_qrrho_checked'
                or result.get('identity') != ident
                or result.get('implementation_certificate_sha256') != proof):
            return None
        execution=result.get('kernel_execution') or {}
        owner=result.get('owner') or {}
        if (execution.get('library_sha256')!=Q.LIBRARY_SHA256
                or execution.get('source_sha256')!=sha(Path(Q.__file__).with_name('xtb_kernel.py'))
                or execution.get('pid')==owner.get('pid') or execution.get('host')!=owner.get('host')
                or set(execution.get('thread_limits',{}).values())!={'1'}
                or not execution.get('cpu_affinity')
                or not set(execution['cpu_affinity'])<=set(owner.get('cpu_affinity',[]))):return None
        resolved={}
        for path, digest in result['artifacts_sha256'].items():
            local=target.with_suffix('')/Path(path).name
            candidate=local if local.is_file() else Path(path)
            if sha(candidate) != digest:return None
            resolved[str(candidate)]=digest
        result=dict(result,artifacts_sha256=resolved,postprocessing_path=str(target))
        from ase.io import read
        atoms = read(geometry)
        matrix = np.load(Path(native_path).with_suffix('') / 'hessian_cartesian.npy')
        recomputed = Q.nuclear_thermals(matrix, atoms.positions, atoms.numbers, result['symmetry_number'])
        for key in ('rigid_rank', 'internal_mode_count', 'retained_soft_mode_indices', 'imaginary_noise_mode_indices'):
            if result.get(key) != recomputed[key]:return None
        if (len(result['frequencies_cm']) != recomputed['internal_mode_count']
                or np.max(abs(np.asarray(result['frequencies_cm']) - recomputed['frequencies_cm']), initial=0) > 1e-6):
            return None
        for cutoff in ('25', '50', '100'):
            value = result['nuclear_by_cutoff_eV'][cutoff]
            if not math.isfinite(value) or abs(value - recomputed['nuclear_by_cutoff_eV'][cutoff]) > 1e-8:
                return None
        comparisons = result['CLI_comparisons']
        if {c['cutoff_cm'] for c in comparisons} != {25, 50, 100}:return None
        serialized_path = next(Path(p) for p in result['artifacts_sha256'] if Path(p).name == 'hessian')
        serialized = _serialized_frequencies(serialized_path, recomputed['atomic_masses_amu'])
        retained = abs(serialized[abs(serialized) > 1])
        rot, average, linear = Q.inertia_properties(atoms.positions, recomputed['atomic_masses_amu'])
        for item in comparisons:
            filtered = Q.kernel_thermal(retained, rot, average, sum(recomputed['atomic_masses_amu']),
                                       result['symmetry_number'], cutoff_cm=item['cutoff_cm'],
                                       linear=linear)['G_Ha'] * Q.HARTREE_EV
            if (item.get('CLI_retained_modes') != len(retained)
                    or item.get('reference') != 'exact serialized Cartesian Hessian and CLI mode filter'
                    or abs(item.get('explicit_filtered_kernel_eV', math.inf) - filtered) > 1e-8
                    or not math.isfinite(item.get('absolute_difference_eV', math.inf))
                    or not 0 <= item['absolute_difference_eV'] <= 1e-6):return None
        nuclear = result['nuclear_by_cutoff_eV']['50']
        spin = -8.617333262145e-5 * 298.15 * math.log(ident['mult'])
        if (abs(result['g_nuclear_qrrho_eV'] - nuclear) > 1e-8
                or abs(result['g_spin_partition_eV'] - spin) > 1e-8
                or abs(result['g_thermal_total_eV'] - nuclear - spin) > 1e-8):return None
        if source_identity(native_path, geometry)[0] != ident:return None
        return result
    except (OSError, ValueError, TypeError, KeyError, RuntimeError, StopIteration):
        return None


def ensure_projected_thermal(native_path, geometry):
    """Process an admitted raw native record without changing its original bytes."""
    native_path, geometry = Path(native_path), Path(geometry)
    proof = certificate()
    target, ident = cache_path(native_path, geometry)
    target.parent.mkdir(exist_ok=True)
    with target.with_suffix('.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        existing = load_projected_thermal(native_path, geometry)
        if existing is not None:return 'present', existing
        if target.exists():
            raise RuntimeError('Invalid projected thermal record must be preserved and inspected')
        work = target.with_suffix('')
        if work.exists():
            raise RuntimeError('Partial projected thermal outputs require ownership review before resuming')
        _, raw, hessian = source_identity(native_path, geometry)
        if (raw.get('status') != 'derivative_checked_solution_thermal'
                or any(raw.get(k) is not True for k in ('converged', 'stationary_minimum', 'grid_response'))
                or raw.get('gradient_max', math.inf) > .00045
                or abs(raw.get('smd_spin_diagnostics', {}).get('spin_contam', math.inf)) > .1):
            raise RuntimeError('Backing native record is not a qualified stationary surface')
        work.mkdir()
        owner = dict(host=os.uname().nodename, pid=os.getpid(),
                     process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                     cpu_affinity=sorted(os.sched_getaffinity(0)), started=time.time())
        atomic_write_json(work / 'owner.json', owner)
        from ase.io import read
        atoms = read(geometry)
        matrix = np.load(hessian)
        masses = Q.atomic_masses(atoms.numbers)
        projected, rank = project_rigid_modes(matrix, atoms.positions, masses)
        np.save(work / 'projected_hessian.npy', projected)
        write_xtb_hessian(work / 'hessian', projected)
        (work / 'mol.xyz').write_bytes(geometry.read_bytes())
        (work / 'qrrho_producer.py').write_bytes(QRRHO_SOURCE_BYTES)
        (work / 'xtb_kernel_producer.py').write_bytes(KERNEL_SOURCE_BYTES)
        (work / 'record_producer.py').write_bytes(PRODUCER_SOURCE_BYTES)
        frequencies, internal_rank = Q.internal_vibrations(matrix, atoms.positions, masses)
        if rank != internal_rank:raise RuntimeError('Rigid-mode classifications disagree')
        rot, average, linear = Q.inertia_properties(atoms.positions, masses)
        serialized = _serialized_frequencies(work / 'hessian', masses)
        comparisons, symmetry_number = [], None
        for cutoff in (25, 50, 100):
            log = work / f'thermo_CLI_{cutoff}.log'
            with log.open('w') as stream:
                process = subprocess.run([str(XTB), 'thermo', '--sthr', str(cutoff),
                                          '--temp', '298.15', 'mol.xyz', 'hessian'], cwd=work,
                                         stdout=stream, stderr=subprocess.STDOUT, timeout=120)
            if process.returncode:raise RuntimeError('Projected CLI reproduction control failed')
            text = log.read_text()
            number = int(re.search(r'rotational number\s+(\d+)', text)[1])
            if symmetry_number is None:symmetry_number = number
            if number != symmetry_number:raise RuntimeError('Symmetry detection changed between cutoffs')
            count = int(re.search(r'# frequencies\s+(\d+)', text)[1])
            retained = abs(serialized[abs(serialized) > 1])
            if count != len(retained):raise RuntimeError('Projected CLI control has unexplained mode count')
            cli = parse_xtb_thermal_eV(text)
            filtered = Q.kernel_thermal(retained, rot, average, sum(masses), number,
                                        cutoff_cm=cutoff, linear=linear)['G_Ha'] * Q.HARTREE_EV
            error = abs(filtered - cli)
            if error > 1e-6:raise RuntimeError('Installed kernel/CLI reproduction failed')
            comparisons.append(dict(cutoff_cm=cutoff, CLI_retained_modes=count,
                                     reference='exact serialized Cartesian Hessian and CLI mode filter',
                                     explicit_filtered_kernel_eV=filtered, CLI_eV=cli,
                                     absolute_difference_eV=error))
        result = Q.nuclear_thermals(matrix, atoms.positions, atoms.numbers, symmetry_number)
        if source_identity(native_path, geometry)[0] != ident:
            raise RuntimeError('Native source changed during projected thermochemistry')
        spin = -8.617333262145e-5 * 298.15 * math.log(ident['mult'])
        result.update(status='projected_internal_qrrho_checked', identity=ident,
                      implementation_certificate_sha256=proof, symmetry_number=symmetry_number,
                      CLI_comparisons=comparisons, g_nuclear_qrrho_eV=result['nuclear_by_cutoff_eV']['50'],
                      g_spin_partition_eV=spin, g_thermal_total_eV=result['nuclear_by_cutoff_eV']['50'] + spin,
                      native_source_path=str(native_path), native_source_hessian=str(hessian),
                      postprocessing_path=str(target), owner=owner, finished=time.time(),
                      kernel_execution=Q.kernel_execution(),
                      producer_hash_convention='source snapshots at module import; kernel source verified by subprocess replies',
                      artifacts_sha256={str(p): sha(p) for p in work.iterdir() if p.is_file()})
        atomic_write_json(target, result)
        loaded = load_projected_thermal(native_path, geometry)
        if loaded is None:raise RuntimeError('Projected partition failed its own independent-input replay')
        return 'computed', loaded


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    selected=parser.add_mutually_exclusive_group(required=True)
    selected.add_argument('--all',action='store_true',help='process all currently admitted backing native Hessians')
    selected.add_argument('--only',help='comma-separated molecule IDs or ID:state selectors')
    parser.add_argument('--include-spin-helpers',action='store_true')
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--out',type=Path,required=True,help='new run directory for ownership and results')
    args=parser.parse_args(argv)
    if args.threads<1:parser.error('positive --threads required')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():
        resources=host_capacity(args.threads)
        os.sched_setaffinity(0,cpu_allocation(args.threads))
    from redox.core.common import DFT,read_manifest
    from redox.qm.native_thermal import _load_native_surface_record,cache_path as native_path
    rows=read_manifest()
    available={token for row in rows for token in (row['id'],row['id']+':'+row['state'])}
    selectors=set(args.only.split(',')) if args.only else available
    if selectors-available:parser.error('unknown declared selectors: '+', '.join(sorted(selectors-available)))
    paths={DFT/row['id']/row['state'] for row in rows
           if row['id'] in selectors or row['id']+':'+row['state'] in selectors}
    if args.include_spin_helpers:
        gids={p.parent.name for p in paths}
        paths.update(p for gid in gids for p in (ROOT/'calcs/spincheck'/gid).glob('*_m*'))
    args.out.mkdir(parents=True,exist_ok=False)
    report=dict(status='running',host=os.uname().nodename,pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                cpu_affinity=sorted(os.sched_getaffinity(0)),resources=resources,
                started=time.time(),results=[],quantum_jobs_submitted=0,raw_records_modified=False)
    destination=args.out/'run.json'
    atomic_write_json(destination,report)
    consecutive=0
    for sd in sorted(paths):
        item=dict(state_dir=str(sd),started=time.time())
        try:
            raw=json.loads((sd/'result.json').read_text()) if (sd/'result.json').is_file() else None
            if raw is None or _load_native_surface_record(sd,raw['charge'],raw['mult']) is None:
                item['status']='skipped_backing_surface_pending'
            else:
                path,_=native_path(sd,raw['charge'],raw['mult'])
                kind,result=ensure_projected_thermal(path,sd/'opt.xyz')
                item.update(status='corrected_partition_checked',kind=kind,
                            partition=result['postprocessing_path'],g_thermal_total_eV=result['g_thermal_total_eV'])
            consecutive=0
        except Exception as exc:
            item.update(status='failed_preserved',error=repr(exc));consecutive+=1
        item['finished']=time.time();report['results'].append(item);atomic_write_json(destination,report)
        print(sd.parent.name,sd.name,item['status'],flush=True)
        if consecutive>=3:
            report['status']='three_failure_circuit_breaker';break
    else:report['status']='bounded_processing_complete'
    report.update(finished=time.time(),checked=sum(x['status']=='corrected_partition_checked' for x in report['results']),
                  skipped=sum(x['status']=='skipped_backing_surface_pending' for x in report['results']),
                  failed=sum(x['status']=='failed_preserved' for x in report['results']))
    atomic_write_json(destination,report)
    if report['failed']:raise SystemExit(1)


if __name__=='__main__':main()
