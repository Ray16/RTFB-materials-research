"""Saved optimizer reuse must fail closed and retain the ordinary energy call."""
import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def resume():
    path = Path(__file__).resolve().parents[1] / 'scripts/validation/resume_benchmark_spin.py'
    spec = importlib.util.spec_from_file_location('spin_resume_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def checkpoint(tmp_path):
    folder = tmp_path / 'checkpoint'; folder.mkdir()
    geometry = '2\ncheckpoint\nH 0 0 0\nH 0 0 0.74\n'
    for name in ('optimizer_returned.xyz', 'last_evaluated.xyz'):
        (folder / name).write_text(geometry)
    (folder / 'optimizer_returned.json').write_text(json.dumps(dict(optimizer_converged=True)))
    (folder / 'latest.json').write_text(json.dumps(dict(SCF_converged=True)))
    return folder


def test_saved_callback_or_changed_return_does_not_establish_convergence(resume, checkpoint):
    (checkpoint / 'optimizer_returned.json').write_text(json.dumps(dict(optimizer_converged=False)))
    with pytest.raises(ValueError, match='real converged'):
        resume.geometry_evidence(checkpoint)
    (checkpoint / 'optimizer_returned.json').write_text(json.dumps(dict(optimizer_converged=True)))
    (checkpoint / 'last_evaluated.xyz').write_text('2\nchanged\nH 0 0 0\nH 0 0 0.9\n')
    with pytest.raises(ValueError, match='does not match'):
        resume.geometry_evidence(checkpoint)


def test_exact_saved_geometry_retains_real_energy_call_and_provenance(resume, checkpoint, tmp_path, monkeypatch):
    from redox.qm import dft
    calls = []
    def energy(xyz, q, m, **kwargs):
        calls.append((Path(xyz).read_bytes(), q, m, kwargs))
        return dict(optimized=False, converged_smd=True, e_smd_eV=-10)
    monkeypatch.setattr(dft, 'dft_smd', energy)
    helper = tmp_path / 'opt.xyz'
    chosen = resume.ROOT / 'calcs/dft' / resume.GID / resume.STATE / 'opt.xyz'
    with resume.reuse_geometry(checkpoint, helper):
        result = dft.dft_smd(chosen, -2, 3, backend='gpu', opt_out=helper, do_freq=False)
    assert calls == [((checkpoint / 'optimizer_returned.xyz').read_bytes(), -2, 3,
                      dict(backend='gpu', do_opt=False, do_freq=False))]
    assert result['optimized'] is True and result['e_smd_eV'] == -10
    assert result['optimization_checkpoint']['SCF_converged'] is True
    assert dft.dft_smd is energy


def test_wrong_spin_or_failed_energy_cannot_install_reused_success(resume, checkpoint, tmp_path, monkeypatch):
    from redox.qm import dft
    def fail(*args, **kwargs):
        raise RuntimeError('SCF did not finish')
    monkeypatch.setattr(dft, 'dft_smd', fail)
    helper = tmp_path / 'opt.xyz'
    chosen = resume.ROOT / 'calcs/dft' / resume.GID / resume.STATE / 'opt.xyz'
    with resume.reuse_geometry(checkpoint, helper), pytest.raises(ValueError, match='declared helper'):
        dft.dft_smd(chosen, -2, 1, backend='gpu', opt_out=helper, do_freq=False)
    with pytest.raises(RuntimeError, match='SCF did not finish'):
        with resume.reuse_geometry(checkpoint, helper):
            dft.dft_smd(chosen, -2, 3, backend='gpu', opt_out=helper, do_freq=False)
    assert dft.dft_smd is fail
    assert not (helper.parent / 'result.json').exists()
