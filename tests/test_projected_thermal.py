"""Addressed partition caches preserve source data and reject stale evidence."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from redox.qm import projected_thermal as P
from redox.qm import qrrho as Q


@pytest.fixture
def source(tmp_path, monkeypatch):
    # Isolate cache semantics from the deployment certificate. The numerical
    # kernel itself is separately tested against xTB's upstream reference.
    monkeypatch.setattr(P, 'certificate', lambda: 'qualified-test-certificate')
    geometry = tmp_path/'opt.xyz'
    geometry.write_text('2\nH2 soft internal-mode fixture\nH -0.5 0 0\nH 0.5 0 0\n')
    native = tmp_path/'native_thermal'/'source.json'
    native.parent.mkdir()
    native.write_text(json.dumps(dict(status='derivative_checked_solution_thermal',
        identity=dict(geometry_sha256=hashlib.sha256(geometry.read_bytes()).hexdigest(), charge=0, mult=1),
        converged=True, stationary_minimum=True, grid_response=True, gradient_max=0.,
        smd_spin_diagnostics=dict(spin_contam=0.))))
    native.with_suffix('').mkdir()
    stretch = np.array([-1., 0., 0., 1., 0., 0.]) / np.sqrt(2)
    mass = Q.atomic_masses([1, 1])[0]
    matrix = (.391 / Q.AUTORCM)**2 * mass * Q.AMUTOAU * np.outer(stretch, stretch)
    np.save(native.with_suffix('')/'hessian_cartesian.npy', matrix)
    return native, geometry


def test_soft_mode_partition_is_cached_without_modifying_source(source):
    native, geometry = source
    files = [native, geometry, native.with_suffix('')/'hessian_cartesian.npy']
    original = {p: p.read_bytes() for p in files}
    kind, result = P.ensure_projected_thermal(native, geometry)
    assert kind == 'computed' and result['internal_mode_count'] == 1
    assert result['retained_soft_mode_indices'] == [0]
    assert result['frequencies_cm'] == pytest.approx([.391])
    assert all(c['CLI_retained_modes'] == 0 for c in result['CLI_comparisons'])
    assert P.ensure_projected_thermal(native, geometry)[0] == 'present'
    assert all(p.read_bytes() == value for p, value in original.items())


@pytest.mark.parametrize('damage', ['total', 'mode_count', 'artifact', 'source', 'kernel_owner'])
def test_stale_or_inconsistent_partition_never_loads(source, damage):
    native, geometry = source
    _, record = P.ensure_projected_thermal(native, geometry)
    target, _ = P.cache_path(native, geometry)
    if damage == 'total':
        record['g_thermal_total_eV'] += .01
        target.write_text(json.dumps(record))
    elif damage == 'mode_count':
        record['internal_mode_count'] += 1
        target.write_text(json.dumps(record))
    elif damage == 'artifact':
        path = next(Path(p) for p in record['artifacts_sha256'] if p.endswith('thermo_CLI_50.log'))
        path.write_text(path.read_text() + '\nchanged evidence\n')
    elif damage == 'kernel_owner':
        record['kernel_execution']['pid']=record['owner']['pid']
        target.write_text(json.dumps(record))
    else:
        geometry.write_text(geometry.read_text().replace('0.5 0 0', '0.6 0 0'))
    assert P.load_projected_thermal(native, geometry) is None
    with pytest.raises(RuntimeError):
        P.ensure_projected_thermal(native, geometry)


def test_missing_implementation_certificate_is_not_an_implicit_admission(source, monkeypatch):
    native, geometry = source
    P.ensure_projected_thermal(native, geometry)
    def missing():raise RuntimeError('qualification missing')
    monkeypatch.setattr(P, 'certificate', missing)
    assert P.load_projected_thermal(native, geometry) is None


def test_atomic_state_promotion_keeps_partition_evidence_readable(source, tmp_path):
    import shutil
    native, geometry = source
    _, result = P.ensure_projected_thermal(native, geometry)
    copied = tmp_path.parent/(tmp_path.name+'_promoted')
    shutil.copytree(tmp_path, copied)
    shutil.rmtree(tmp_path)
    new_native, new_geometry = copied/native.relative_to(tmp_path), copied/'opt.xyz'
    loaded = P.load_projected_thermal(new_native, new_geometry)
    assert loaded and loaded['g_thermal_total_eV']==result['g_thermal_total_eV']
    assert loaded['postprocessing_path'].startswith(str(copied))
    assert all(Path(p).is_file() for p in loaded['artifacts_sha256'])
