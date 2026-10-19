"""Source-bound equilibrium free energies for two qualified spin surfaces.

Component native free energies already include spin degeneracy. This module never
runs quantum calculations and does not establish spin-conversion kinetics.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

MODEL = "two-spin equilibrium; component spin degeneracy already included"
TEMPERATURE_K = 298.15
POLICY_NAME = "spin_ensemble_policy.json"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def equilibrium(chosen_G, alternative_G, temperature_K=TEMPERATURE_K):
    if (any(type(x) not in (int, float) or not math.isfinite(x)
            for x in (chosen_G, alternative_G, temperature_K)) or temperature_K <= 0):
        raise ValueError("Finite component free energies and positive temperature required")
    kt = 8.617333262145e-5 * temperature_K
    gap = alternative_G - chosen_G
    effective = min(chosen_G, alternative_G) - kt * math.log1p(math.exp(-abs(gap)/kt))
    population = math.exp(-max(gap, 0)/kt) / (1 + math.exp(-abs(gap)/kt))
    return dict(free_energy_eV=effective, correction_eV=effective-chosen_G,
                alternative_population=population, gap_eV=gap)


def _current_components(state_dir, helper_dir):
    """Read actual active energies, qualified partitions and exact molecular identities."""
    from redox.core.common import ROOT, DFT, smiles_to_xyz_map
    from redox.core.protocol import load_record, record_path
    from redox.qm.dft import _opt_protocol_mismatch, _uma_mult
    from redox.qm.native_thermal import load_native_thermal, cache_path, NATIVE_THERMAL
    from redox.validation.spins import declared, snapshot, alternate_mult

    state_dir, helper_dir = Path(state_dir).resolve(), Path(helper_dir).resolve()
    gid, state = state_dir.parent.name, state_dir.name
    if state_dir != (DFT/gid/state).resolve():
        raise ValueError("Ensemble selection requires a canonical declared state")
    identity = snapshot(gid, state)
    if identity is None:
        raise ValueError("Selected active energy is missing")
    q, mult = identity['charge'], identity['mult']
    alt = alternate_mult(mult)
    if helper_dir != (ROOT/'calcs/spincheck'/gid/f'{state}_m{alt}').resolve():
        raise ValueError("Alternative helper identity differs from the declared spin pair")
    row = declared(gid, state)
    if int(row['charge']) != q or _uma_mult(gid,state,int(row.get('mult_hint') or row.get('mult') or 1)) != mult:
        raise ValueError("Selected charge or multiplicity differs from its declared decision")
    smiles = row['smiles']
    components, evidence, thermal_protocol, seed_identity = [], {}, None, None
    for role, directory, m in [('chosen', state_dir, mult), ('alternative', helper_dir, alt)]:
        raw = json.loads((directory/'result.json').read_text())
        if (type(raw.get('charge')) is not int or type(raw.get('mult')) is not int
                or raw.get('charge') != q or raw.get('mult') != m or raw.get('optimized') is not True
                or _opt_protocol_mismatch(raw, q)
                or smiles_to_xyz_map(smiles, directory/'opt.xyz') is None):
            raise ValueError("Spin component geometry, optimization protocol or exact graph failed")
        if role == 'alternative':
            seed_identity = raw.get('spin_source_identity')
            # This is a new thermochemical comparison of two independently
            # qualified own minima, not reuse of an old electronic-gap report.
            # The starting geometry may differ after a canonical repair; its
            # chemistry, selected charge/spin and Hamiltonian must still match.
            if (not isinstance(seed_identity, dict) or any(seed_identity.get(k) != identity[k]
                    for k in ('id','state','charge','mult','smiles_sha256','protocol'))):
                raise ValueError("Alternative helper belongs to another chemical surface")
        policy_path = directory/'thermal_policy.json'
        selected_thermal = json.loads(policy_path.read_text())
        if (selected_thermal.get('thermal_definition') != 'native_solution_dft_qrrho_spin'
                or selected_thermal.get('native_protocol') != NATIVE_THERMAL):
            raise ValueError("Both components must select the qualified native thermal model")
        energy = load_record(directory, q, m)
        native = load_native_thermal(directory, q, m)
        if energy is None or native is None:
            raise ValueError("Both spin components require active energies and corrected native partitions")
        contamination = energy.get('smd_spin_contam')
        if (type(contamination) not in (int, float) or not math.isfinite(contamination)
                or abs(contamination) > .1):
            raise ValueError("Spin component contamination is unqualified")
        protocol = native['thermal_partition_protocol']
        if thermal_protocol is not None and protocol != thermal_protocol:
            raise ValueError("Spin components have different nuclear partition protocols")
        if protocol.get('temperature_K') != TEMPERATURE_K:
            raise ValueError("Native component temperature differs from the ensemble model")
        thermal_protocol = protocol
        raw_native = cache_path(directory, q, m)[0]
        partition = Path(native['thermal_postprocessing_path'])
        required = [directory/'result.json', directory/'opt.xyz', policy_path, record_path(directory, q, m),
                    raw_native, raw_native.with_suffix('')/'hessian_cartesian.npy', partition]
        # The corrected reader has checked artifact identities; bind the same local
        # artifacts so later edits cannot retain a valid ensemble selection.
        post = json.loads(partition.read_text())
        for path, expected in post['artifacts_sha256'].items():
            original = Path(path)
            local = partition.with_suffix('')/original.name
            actual = local if local.is_file() else original
            if digest(actual) != expected:
                raise ValueError("Corrected partition artifact changed")
            required.append(actual)
        evidence.update({str(p): digest(p) for p in required})
        e, g = energy['e_smd_eV'], native['g_thermal_total_eV']
        components.append(dict(role=role, directory=str(directory), charge=q, mult=m,
                               electronic_eV=e, thermal_eV=g, free_energy_eV=e+g,
                               geometry_sha256=digest(directory/'opt.xyz'),
                               energy_path=str(record_path(directory,q,m)),
                               native_path=str(raw_native), partition_path=str(partition)))
    if components[1]['electronic_eV'] < components[0]['electronic_eV'] - 1e-4:
        raise ValueError("A lower electronic spin must complete canonical promotion first")
    return dict(chosen_identity=identity, alternative_seed_identity=seed_identity,
                helper=str(helper_dir), components=components,
                native_partition_protocol=thermal_protocol, evidence_sha256=evidence)


def prepare_policy(state_dir, helper_dir):
    """Build an unselected policy from existing evidence, with no canonical writes."""
    current = _current_components(state_dir, helper_dir)
    a, b = current['components']
    mixed = equilibrium(a['free_energy_eV'], b['free_energy_eV'])
    return dict(schema=1, model=MODEL, temperature_K=TEMPERATURE_K,
                **current, equilibrium=mixed, single_surface_kinetics_validated=False)


def evaluate_policy(state_dir, policy):
    """Recompute and validate a prepared/selected policy; never trust saved populations."""
    if (type(policy.get('schema')) is not int or policy.get('schema') != 1 or policy.get('model') != MODEL
            or policy.get('temperature_K') != TEMPERATURE_K
            or policy.get('single_surface_kinetics_validated') is not False):
        raise ValueError("Unsupported spin ensemble model or kinetic claim")
    current = _current_components(state_dir, policy['helper'])
    for key in ('chosen_identity', 'alternative_seed_identity', 'helper', 'components',
                'native_partition_protocol', 'evidence_sha256'):
        if policy.get(key) != current[key]:
            raise ValueError(f"Spin ensemble dependency changed: {key}")
    chosen, alternative = current['components']
    expected = equilibrium(chosen['free_energy_eV'], alternative['free_energy_eV'])
    saved = policy.get('equilibrium') or {}
    for key, value in expected.items():
        actual = saved.get(key)
        if type(actual) not in (int, float) or not math.isfinite(actual) or abs(actual-value) > 1e-8:
            raise ValueError("Stored ensemble arithmetic differs from current components")
    if any(digest(p) != h for p, h in current['evidence_sha256'].items()):
        raise ValueError("Spin ensemble input changed during evaluation")
    return dict(expected, helper=current['helper'], components=current['components'],
                chosen_identity=current['chosen_identity'],
                evidence_sha256=current['evidence_sha256'])


def apply_selected_ensemble(record, state_dir):
    """Overlay a selected equilibrium correction, or withhold a stale selection."""
    record = dict(record)
    # Producer JSON flags cannot pretend that a selection was verified by this reader.
    for key in list(record):
        if key.startswith('spin_ensemble_') or key == 'single_surface_kinetics_validated':
            record.pop(key)
    record['spin_ensemble_status'] = 'not_selected'
    path = Path(state_dir)/POLICY_NAME
    if not path.exists():
        return record
    record['spin_ensemble_policy_path'] = str(path)
    if record.get('thermal_qc') != 'ok' or record.get('native_thermal_record') is None:
        record.update(g_thermal_eV=None, thermal_qc='spin_ensemble_component_unqualified',
                      spin_ensemble_status='invalid')
        return record
    try:
        policy = json.loads(path.read_text())
        checked = evaluate_policy(state_dir, policy)
        selected = policy['components'][0]
        g = record['native_thermal_record']['g_thermal_total_eV']
        if (abs(g-selected['thermal_eV']) > 1e-8
                or abs(record['e_smd_eV']-selected['electronic_eV']) > 1e-8):
            raise ValueError("Primary state differs from selected ensemble component")
        receipt = digest(path)
        record.update(g_thermal_eV=g+checked['correction_eV'],
                      thermal_method='solution DFT two-spin equilibrium; native spin degeneracy included',
                      thermal_protocol_hash='ensemble:'+receipt,
                      spin_ensemble_status='applied', spin_ensemble_policy_sha256=receipt,
                      spin_ensemble_correction_eV=checked['correction_eV'],
                      spin_ensemble_free_energy_eV=checked['free_energy_eV'],
                      spin_ensemble_selected_surface_G_eV=selected['free_energy_eV'],
                      spin_ensemble_helper=checked['helper'],
                      spin_ensemble_alternative_population=checked['alternative_population'],
                      single_surface_kinetics_validated=False)
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        record.update(g_thermal_eV=None, thermal_qc='stale_or_invalid_spin_ensemble',
                      spin_ensemble_status='invalid')
    return record


def admitted_ensemble(record, state_dir):
    """Confirm primary equilibrium use before overriding legacy spin-report gates."""
    if not record or record.get('spin_ensemble_status') != 'applied':
        return None
    path=Path(state_dir)/POLICY_NAME
    try:
        if (not path.is_file() or digest(path) != record.get('spin_ensemble_policy_sha256')
                or record.get('single_surface_kinetics_validated') is not False):
            return None
        checked=evaluate_policy(state_dir,json.loads(path.read_text()))
        if record.get('spin_ensemble_helper')!=checked['helper']:
            return None
        values=[record.get(k) for k in ('e_smd_eV','g_thermal_eV','spin_ensemble_alternative_population')]
        if any(type(v) not in (int,float) or not math.isfinite(v) for v in values):
            return None
        if abs(record['e_smd_eV']+record['g_thermal_eV']-checked['free_energy_eV'])>1e-8:
            return None
        if abs(record['spin_ensemble_alternative_population']-checked['alternative_population'])>1e-8:
            return None
        return checked
    except (OSError,ValueError,KeyError,TypeError,AttributeError):
        return None


def _independent_thermal_inputs(state_dir):
    import importlib.util
    from redox.core.common import ROOT
    spec=importlib.util.spec_from_file_location('selected_equilibrium_independent',
            ROOT/'scripts/validation/audit/recompute_axes.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    raw=json.loads((Path(state_dir)/'result.json').read_text())
    return audit.thermal_inputs(Path(state_dir),raw)


def select_policy(state_dir, policy):
    """Select a prepared policy only after both real readers admit its free energy.

    No quantum calculation runs here. Old selections are retained by hash and a
    failed validation restores the previous selection before releasing writer locks.
    """
    import fcntl
    import os
    from contextlib import ExitStack
    from redox.core import common
    from redox.core.io import atomic_write_json

    if not hasattr(common,'_apply_single_surface_thermal_policy'):
        raise RuntimeError('Shared equilibrium reader hook is not activated')
    directory=Path(state_dir).resolve();gid,state=directory.parent.name,directory.name
    lock_root=common.ROOT/'logs/pipeline_locks';lock_root.mkdir(parents=True,exist_ok=True)
    with ExitStack() as stack:
        for suffix in ('','_native_thermal','_spincheck'):
            lock=stack.enter_context((lock_root/f'{gid}{suffix}.lock').open('a+'))
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        checked=evaluate_policy(directory,policy)
        target=directory/POLICY_NAME
        previous=target.read_bytes() if target.exists() else None
        history=directory/'spin_ensemble_history';history.mkdir(exist_ok=True)
        if previous is not None:
            old=history/(hashlib.sha256(previous).hexdigest()+'.json')
            if old.exists() and old.read_bytes()!=previous:raise RuntimeError('Archived selection collision')
            old.write_bytes(previous)
        atomic_write_json(target,policy)
        new_bytes=target.read_bytes()
        try:
            record=common.read_result(gid,state)
            admitted=admitted_ensemble(record,directory)
            if admitted is None:raise RuntimeError('Primary reader did not apply the selected equilibrium policy')
            thermal,imaginary,definition=_independent_thermal_inputs(directory)
            if (thermal is None or imaginary!=0 or definition!='native_solution_dft_qrrho_spin'
                    or abs(record['e_smd_eV']+thermal-checked['free_energy_eV'])>1e-8):
                raise RuntimeError('Independent reader did not admit the selected equilibrium value')
            if any(digest(p)!=h for p,h in checked['evidence_sha256'].items()):
                raise RuntimeError('Component data changed while selecting equilibrium policy')
        except BaseException:
            failed=history/(hashlib.sha256(new_bytes).hexdigest()+'.failed.json')
            failed.write_bytes(new_bytes)
            if previous is None:target.unlink()
            else:
                restored=target.with_name(f'.{target.name}.{os.getpid()}.restore.tmp')
                restored.write_bytes(previous)
                restored.chmod(target.stat().st_mode & 0o777)
                restored.replace(target)
            raise
        return dict(status='equilibrium_policy_selected_and_both_readers_verified',
                    path=str(target),sha256=digest(target),free_energy_eV=checked['free_energy_eV'],
                    alternative_population=checked['alternative_population'],single_surface_kinetics_validated=False)
