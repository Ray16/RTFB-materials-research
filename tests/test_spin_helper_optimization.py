"""A stopped helper optimizer must not publish a replacement geometry."""
import json
from types import SimpleNamespace

import pytest

from redox.qm import dft
from redox.validation import spins


@pytest.mark.parametrize('geometry_ok,scf_ok', [(True,True),(False,True),(True,False)])
def test_helper_retains_optimizer_and_scf_convergence(monkeypatch,geometry_ok,scf_ok):
    scanner=SimpleNamespace(converged=scf_ok)
    gradient=SimpleNamespace(grid_response=False,as_scanner=lambda:scanner)
    mf=SimpleNamespace(nuc_grad_method=lambda:gradient)
    final_mol=object()

    def kernel(actual,**kwargs):
        assert actual is scanner and gradient.grid_response is True
        assert kwargs['assert_convergence'] is True
        assert kwargs['maxsteps']==120
        assert kwargs['convergence_gmax']==5e-5
        return geometry_ok,final_mol

    monkeypatch.setattr(dft,'_geomopt_fn',lambda backend:kernel)
    if not (geometry_ok and scf_ok):
        with pytest.raises(RuntimeError,match='did not converge'):
            spins._optimize_helper_minimum(mf)
    else:
        mol,proof=spins._optimize_helper_minimum(mf)
        assert mol is final_mol
        assert proof==dict(optimized=True,geometry_optimization_converged=True,
                           optimization_SCF_converged=True,opt_gradient_grid_response=True)


@pytest.mark.parametrize('geometry_ok,scf_ok',[(False,True),(True,False)])
def test_failed_helper_refinement_preserves_input_and_writes_no_geometry(
    tmp_path,monkeypatch,geometry_ok,scf_ok
):
    helper=tmp_path/'helper';helper.mkdir()
    xyz=helper/'opt.xyz';xyz.write_text('2\noriginal\nH 0 0 0\nH 0 0 .74\n')
    result=helper/'result.json';result.write_text(json.dumps(dict(charge=0,mult=1)))
    original={p.name:p.read_bytes() for p in helper.iterdir()}
    monkeypatch.setattr(spins,'BASE',tmp_path/'review')
    def unqualified(*args,**kwargs):
        raise RuntimeError('Source is not stationary')
    monkeypatch.setattr(spins,'ensure_native_thermal',unqualified)
    gradient=SimpleNamespace(grid_response=False,as_scanner=lambda:SimpleNamespace(converged=scf_ok))
    mf=SimpleNamespace(converged=True,kernel=lambda:None,nuc_grad_method=lambda:gradient,
                       grids=SimpleNamespace(level=3))
    monkeypatch.setattr(dft,'_build_smd_mf',lambda *a,**k:mf)
    monkeypatch.setattr(dft,'_geomopt_fn',lambda backend:lambda *a,**k:(geometry_ok,None))
    with pytest.raises(RuntimeError,match='Bounded helper minimum repair failed'):
        spins._helper_native('test_parent','neu',helper)
    assert original=={p.name:p.read_bytes() for p in helper.iterdir()}
    assert not list((tmp_path/'review').glob('**/branch_*/opt.xyz'))
    reports=list((tmp_path/'review').glob('**/status.json'))
    assert len(reports)==1
    assert 'did not converge' in json.loads(reports[0].read_text())['branches'][0]['error']
