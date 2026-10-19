"""Solver failures and legacy files cannot masquerade as qualified gas references."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from redox.core.common import ROOT
from redox.qm import gas_reference as gas


@pytest.mark.parametrize('converged,value',[(False,-1.),(True,float('nan')),
    (True,float('inf')),(True,-float('inf')),(True,True),(True,np.bool_(True)),(True,'-1')])
def test_bad_scf_energies_are_refused(converged,value):
    with pytest.raises(RuntimeError):
        gas.energy_Ha(SimpleNamespace(kernel=lambda:value,converged=converged))


def test_energy_remains_in_hartree():
    assert gas.energy_Ha(SimpleNamespace(kernel=lambda:np.float64(-1.25),converged=True))==-1.25


@pytest.mark.parametrize('geometry_ok,scf_ok,finite',[(True,True,True),(False,True,True),
                                                   (True,False,True),(True,True,False)])
def test_gas_optimizer_checks_both_convergence_flags_and_coordinates(monkeypatch,geometry_ok,scf_ok,finite):
    from pyscf.geomopt import geometric_solver
    scanner=SimpleNamespace(converged=scf_ok)
    gradient=SimpleNamespace(grid_response=False,as_scanner=lambda:scanner)
    mf=SimpleNamespace(nuc_grad_method=lambda:gradient)
    coordinates=np.zeros((2,3));coordinates[0,0]=0 if finite else np.nan
    mol=SimpleNamespace(atom_coords=lambda:coordinates)
    def kernel(actual,**kwargs):
        assert actual is scanner and gradient.grid_response is True
        assert kwargs==dict(maxsteps=7,assert_convergence=True)
        return geometry_ok,mol
    monkeypatch.setattr(geometric_solver,'kernel',kernel)
    if geometry_ok and scf_ok and finite:
        assert gas.optimized_molecule(mf,maxsteps=7) is mol
    else:
        with pytest.raises(RuntimeError):gas.optimized_molecule(mf,maxsteps=7)


def test_gas_helper_does_not_silently_accept_a_solution_method():
    with pytest.raises(ValueError,match='gas reference'):
        gas.optimized_molecule(SimpleNamespace(with_solvent=object()))


class DeviceArray:
    def __init__(self,values):self.values=np.asarray(values)
    def __array__(self,*args,**kwargs):raise TypeError('implicit device conversion prohibited')
    def get(self):return self.values.copy()


@pytest.mark.parametrize('unrestricted',[False,True])
def test_homo_transfers_device_arrays_explicitly(unrestricted):
    if unrestricted:
        energy=(DeviceArray([-.3,-.1,.2]),DeviceArray([-.4,-.2,.3]))
        occ=(DeviceArray([1,1,0]),DeviceArray([1,0,0]))
    else:
        energy=DeviceArray([-.3,-.1,.2]);occ=DeviceArray([2,2,0])
    assert gas.homo_eV(SimpleNamespace(mo_energy=energy,mo_occ=occ))==pytest.approx(-.1*gas.HARTREE_EV)


@pytest.mark.parametrize('energy,occ',[([np.nan,.2],[1,0]),([-.1,.2],[0,0]),
                                     ([-.1,.2],[np.nan,0]),([-.1,.2],[1])])
def test_unavailable_or_invalid_homo_is_unknown(energy,occ):
    assert gas.homo_eV(SimpleNamespace(mo_energy=energy,mo_occ=occ)) is None


@pytest.mark.parametrize('corruption',['legacy','method','geometry_input','producer','status','malformed'])
def test_cache_is_qualified_or_preserved_without_recomputation(tmp_path,corruption):
    producer=tmp_path/'producer.py';producer.write_text('original producer')
    identity=dict(smiles='CC',basis='sto-3g',phase='gas')
    path=tmp_path/'case.json'
    record=dict(status='ok',solver_verification=gas.verification(identity,producer))
    if corruption=='legacy':record.pop('solver_verification')
    elif corruption=='method':record['solver_verification']['identity']['basis']='different'
    elif corruption=='geometry_input':record['solver_verification']['identity']['smiles']='CO'
    elif corruption=='producer':producer.write_text('changed producer')
    elif corruption=='status':record['status']='partial'
    path.write_text('[]' if corruption=='malformed' else json.dumps(record))
    before=path.read_bytes()
    with pytest.raises(RuntimeError,match='preserved'):gas.existing_result(path,identity,producer)
    assert path.read_bytes()==before


def test_exact_solver_cache_is_reused_and_duplicate_writers_are_excluded(tmp_path):
    producer=tmp_path/'producer.py';producer.write_text('producer')
    path=tmp_path/'case.json';identity=dict(smiles='CC')
    assert gas.existing_result(path,identity,producer) is None
    record=dict(status='ok',solver_verification=gas.verification(identity,producer))
    path.write_text(json.dumps(record))
    assert gas.existing_result(path,identity,producer)==record
    with gas.record_guard(path):
        with pytest.raises(RuntimeError,match='live writer'):
            with gas.record_guard(path):pytest.fail('second writer admitted')


@pytest.mark.parametrize('value',[None,float('nan'),float('inf'),True,'0.5'])
def test_cache_requires_finite_calculated_outputs(tmp_path,value):
    producer=tmp_path/'producer.py';producer.write_text('producer')
    path=tmp_path/'case.json';identity=dict(smiles='CC')
    path.write_text(json.dumps(dict(status='ok',our_lambda=value,
                                   solver_verification=gas.verification(identity,producer))))
    before=path.read_bytes()
    with pytest.raises(RuntimeError,match='Incomplete reference output preserved'):
        gas.existing_result(path,identity,producer,required_fields=('our_lambda',))
    assert path.read_bytes()==before


WORKERS=[ROOT/'scripts/validation/reorg_d3tales/validate_reorg_worker_d3tales.py',
         ROOT/'scripts/validation/reorg_literature/validate_reorg_worker_inner_acene.py',
         ROOT/'scripts/validation/reorg_literature/validate_reorg_worker_generic.py']


@pytest.fixture(params=WORKERS,ids=['d3tales','acene','generic'])
def worker(request,tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('gas_reference_worker_test',request.param)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    monkeypatch.setattr(module,'CALC',tmp_path/'calc')
    if hasattr(module,'GEO'):monkeypatch.setattr(module,'GEO',tmp_path/'geom')
    return module


def row():
    return dict(id='case',smiles='CC',omega=.3,lit_hole_eV=.5,lit_eV=.5,
                d3_hole=.5,d3_electron=np.nan,n_atoms=8)


def run_case(worker,source,couple='hole'):
    return worker.process(source,couple,'cpu') if 'generic' in worker.__file__ else worker.process(source,'cpu')


def test_each_worker_uses_strict_scf_admission(worker,monkeypatch):
    monkeypatch.setattr(worker,'_mol',lambda *a:None)
    monkeypatch.setattr(worker,'_mf',lambda *a:SimpleNamespace(kernel=lambda:-1.,converged=False))
    args=([],0,0,.3,'cpu') if 'd3tales.py' in worker.__file__ else ([],0,0,'cpu')
    with pytest.raises(RuntimeError,match='unconverged'):worker._E(*args)


def test_each_worker_preserves_legacy_records_before_any_calculation(worker,monkeypatch):
    worker.CALC.mkdir();p=worker.CALC/'case.json';p.write_text('{"status":"ok","our_hole":0.5}')
    before=p.read_bytes()
    monkeypatch.setattr(worker,'_embed',lambda *a:pytest.fail('legacy case recomputed'))
    with pytest.raises(RuntimeError,match='preserved'):run_case(worker,row())
    assert p.read_bytes()==before


def setup_fake_success(worker,monkeypatch):
    monkeypatch.setattr(worker,'_embed',lambda *a:[])
    monkeypatch.setattr(worker,'_opt',lambda *a:[])
    monkeypatch.setattr(worker,'_E',lambda *a,**k:-1.)
    if hasattr(worker,'_couple'):monkeypatch.setattr(worker,'_couple',lambda *a,**k:(.5,{}) if k.get('return_evidence') else .5)
    if hasattr(worker,'_hole_lambda'):monkeypatch.setattr(worker,'_hole_lambda',lambda *a,**k:(.5,.2,.3,{}) if k.get('return_evidence') else (.5,.2,.3))
    if hasattr(worker,'_four_point'):monkeypatch.setattr(worker,'_four_point',lambda *a,**k:(.5,.2,.3,-.1,{}) if k.get('return_evidence') else (.5,.2,.3,-.1))


def test_each_worker_preserves_raw_points_and_correct_four_point_units(worker,monkeypatch):
    neutral=[('H',(0.,0.,0.)),('H',(0.,0.,.7))]
    ion=[('H',(0.,0.,0.)),('H',(0.,0.,1.))]
    monkeypatch.setattr(worker,'_opt',lambda *a:ion)
    def point(atoms,q,spin,*args,**kwargs):
        energy={('neutral',0):-10.,('ion',0):-9.98,
                ('neutral',1):-10.06,('ion',1):-10.1}['neutral' if atoms is neutral else 'ion',q]
        return (energy,-.1) if kwargs.get('want_homo') else energy
    monkeypatch.setattr(worker,'_E',point)
    if hasattr(worker,'_couple'):
        values=worker._couple(neutral,-10.,1,1,.3,'cpu',return_evidence=True)
    elif hasattr(worker,'_hole_lambda'):
        values=worker._hole_lambda(neutral,-10.,'cpu',return_evidence=True)
    else:
        values=worker._four_point(neutral,-10.,1,1,'cpu',return_evidence=True)
    assert values[0]==pytest.approx(1.6327,abs=1e-8)
    evidence=values[-1]
    assert evidence['lambda_unrounded_eV']==pytest.approx(.06*gas.HARTREE_EV)
    assert evidence['neutral_geometry_A']==neutral and evidence['ion_geometry_A']==ion
    assert evidence['energies_Ha']==dict(neutral_at_neutral=-10.,ion_at_ion=-10.1,
                                        ion_at_neutral=-10.06,neutral_at_ion=-9.98)


def test_each_worker_records_actual_method_and_only_reuses_qualified_results(worker,monkeypatch):
    setup_fake_success(worker,monkeypatch)
    worker.BASIS='6-311g(d,p)'
    assert run_case(worker,row())=='ok'
    record=json.loads((worker.CALC/'case.json').read_text())
    assert record['solver_verification']['identity']['basis']=='6-311g(d,p)'
    if 'inner_acene' in worker.__file__:assert '6-311g(d,p)' in record['protocol']
    monkeypatch.setattr(worker,'_embed',lambda *a:pytest.fail('qualified cache recomputed'))
    assert run_case(worker,row())=='skip'


@pytest.mark.parametrize('homo,expected',[(None,'binding_unverified'),(.1,'anion_unbound')])
def test_electron_binding_cannot_be_inferred_from_missing_or_positive_homo(tmp_path,monkeypatch,homo,expected):
    spec=importlib.util.spec_from_file_location('gas_binding_test',WORKERS[2]);worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    monkeypatch.setattr(worker,'CALC',tmp_path/'calc');setup_fake_success(worker,monkeypatch)
    monkeypatch.setattr(worker,'_four_point',lambda *a,**k:(.5,.2,.3,homo,{}) if k.get('return_evidence') else (.5,.2,.3,homo))
    assert run_case(worker,row(),couple='electron')==expected
    record=json.loads((worker.CALC/'case.json').read_text())
    assert record['our_lambda'] is None and record['diagnostic_lambda_eV']==.5
    assert record['anion_unbound'] is (None if homo is None else True)
    assert 'solver_verification' not in record


def test_each_worker_cli_returns_failure_for_incomplete_outputs(worker,tmp_path,monkeypatch,capsys):
    import pandas as pd
    import sys
    monkeypatch.setattr(worker,'ROOT',tmp_path)
    monkeypatch.setattr(worker.pd,'read_csv',lambda *a:pd.DataFrame([row()]))
    def failed(source,*args):
        worker.CALC.mkdir(parents=True,exist_ok=True)
        (worker.CALC/'case.json').write_text('{"status":"error"}')
        return 'error'
    monkeypatch.setattr(worker,'process',failed)
    argv=['worker','--source','unused','--backend','cpu']
    if 'generic' in worker.__file__:argv+=['--tag','test','--couple','hole']
    elif 'inner_acene' in worker.__file__:argv+=['--tag','test']
    else:argv+=['--out',str(tmp_path/'out')]
    monkeypatch.setattr(sys,'argv',argv)
    with pytest.raises(SystemExit,match='failed or remain unqualified'):worker.main()
    assert 'DONE_MARKER' not in capsys.readouterr().out
