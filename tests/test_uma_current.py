"""Default engine routing, ASE state-cache isolation and graph-bound warm starts."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes

from redox.qm import uma_current as current


class Fixed(Calculator):
    implemented_properties=['energy','forces']
    def __init__(self,charge,mult):
        super().__init__();self.charge,self.mult=charge,mult
    def calculate(self,atoms=None,properties=('energy','forces'),system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        self.results=dict(energy=float(10*self.charge+self.mult),forces=np.zeros((len(atoms),3)))


def test_identical_coordinates_with_changed_charge_or_spin_do_not_reuse_cache():
    calls=[]
    def factory(numbers,charge,mult):
        calls.append((numbers,charge,mult));return Fixed(charge,mult)
    calc=current.state_aware_calculator(factory,'legacy')
    atoms=Atoms('OH2',positions=[[0,0,0],[0,.8,.6],[0,-.8,.6]])
    atoms.calc=calc
    for charge,mult,energy in [(0,1,1),(1,2,12),(0,3,3),(0,1,1)]:
        atoms.info.update(charge=charge,spin=mult)
        assert atoms.get_potential_energy()==energy
        np.testing.assert_array_equal(atoms.get_forces(),np.zeros((3,3)))
    assert len(calls)==calc.evaluations==calc.transitions==4
    assert atoms.get_potential_energy()==1 and calc.evaluations==4


def test_spin_scan_resolves_nonzero_gap_at_identical_stationary_geometries(tmp_path):
    from ase.io import write
    from redox.qm.uma import relax_state_scan
    xyz=tmp_path/'water.xyz';write(xyz,Atoms('OH2',positions=[[0,0,0],[0,.8,.6],[0,-.8,.6]]))
    calc=current.state_aware_calculator(lambda numbers,q,m:Fixed(q,m),'legacy')
    result=relax_state_scan(calc,[xyz],0,1,.01,5)
    assert result['chosen_mult']==1 and result['mults_scanned']==[1,3]
    assert result['spin_gap_eV']==2 and result['converged']


def test_invalid_explicit_state_stops_before_a_predictor_is_created():
    def forbidden(*args):raise AssertionError('Model must not load')
    atoms=Atoms('OH2');atoms.info.update(charge=0,spin=2)
    atoms.calc=current.state_aware_calculator(forbidden,'legacy')
    with pytest.raises(ValueError):atoms.get_potential_energy()


def test_explicit_cache_reset_recomputes_the_underlying_prediction():
    atoms=Atoms('OH2');atoms.info.update(charge=0,spin=1)
    calc=current.state_aware_calculator(lambda numbers,q,m:Fixed(q,m),'legacy')
    atoms.calc=calc
    assert atoms.get_potential_energy()==1
    assert calc.fixed.results
    calc.reset()
    assert not calc.fixed.results
    assert atoms.get_potential_energy()==1 and calc.evaluations==2


def test_evaluation_checkpoint_keeps_explicit_state_and_unqualified_coordinates(tmp_path):
    atoms=Atoms('OH2');atoms.info.update(charge=0,spin=1)
    atoms.calc=current.state_aware_calculator(lambda numbers,q,m:Fixed(q,m),'legacy',tmp_path)
    atoms.get_potential_energy()
    record=json.loads((tmp_path/'state_0001/latest.json').read_text())
    assert record['charge']==0 and record['mult']==1 and record['numbers']==[8,1,1]
    assert record['optimized_or_qualified'] is False
    assert record['coordinates_A']==atoms.positions.tolist()


def test_false_nitro_fallback_is_rejected_after_evaluation():
    def factory(numbers,q,m):
        calc=Fixed(q,m)
        calc.predictor=SimpleNamespace(inference_settings=SimpleNamespace(nitro=False))
        return calc
    atoms=Atoms('OH2');atoms.info.update(charge=0,spin=1)
    atoms.calc=current.state_aware_calculator(factory,'nitro')
    with pytest.raises(RuntimeError,match='fell back'):atoms.get_potential_energy()


@pytest.mark.parametrize('capability',[(8,0),(8,9),(9,0),(10,0),(12,0)])
def test_nitro_capable_architecture_selects_nitro(capability):
    assert current.select_backend(capability)=='nitro'


@pytest.mark.parametrize('capability',[(7,0),(7,5),(6,1)])
def test_stock_is_the_default_on_every_non_nitro_capable_card(capability):
    # Stock by default: V100 (7.0), T4 (7.5) and anything below cc 8 use legacy.
    assert current.select_backend(capability)=='legacy'


def test_require_nitro_cannot_fall_back_on_a_non_capable_card():
    with pytest.raises(RuntimeError):current.select_backend((7,0),require_nitro=True)
    with pytest.raises(RuntimeError):current.select_backend((7,5),require_nitro=True)


def test_partial_preferred_cache_never_silently_selects_legacy(tmp_path,monkeypatch):
    monkeypatch.setattr(current,'ROOT',tmp_path)
    policy=dict(nitro_output='nitro')
    assert current.nitro_state_directory('g','s',policy) is None
    (tmp_path/'nitro/g/s').mkdir(parents=True)
    with pytest.raises(RuntimeError,match='partial'):current.nitro_state_directory('g','s',policy)


def test_valid_cache_is_bound_to_weights_references_geometry_and_manifest(tmp_path,monkeypatch):
    from ase.io import write
    from redox.core import common
    monkeypatch.setattr(current,'ROOT',tmp_path)
    policy_path=tmp_path/'policy.json';policy_path.write_text('{}')
    monkeypatch.setattr(current,'POLICY',policy_path)
    policy=dict(nitro_output='nitro',checkpoint_sha256='weights',reference_sha256={'atom_refs':'refs'})
    folder=tmp_path/'nitro/g/s';folder.mkdir(parents=True)
    write(folder/'relaxed.xyz',Atoms('OH2',positions=[[0,0,0],[0,.8,.6],[0,-.8,.6]]))
    record=dict(converged=True,id='g',state='s',smiles='O',charge=0,chosen_mult=1,
        geometry_sha256=current.digest(folder/'relaxed.xyz'),
        model_identity=dict(checkpoint_sha256='weights',reference_sha256={'atom_refs':'refs'}),
        UMA_execution=dict(actual_mode='nitro',policy_sha256=current.digest(policy_path)))
    monkeypatch.setattr(common,'read_manifest',lambda:[dict(id='g',state='s',smiles='O',charge='0')])
    (folder/'result.json').write_text(json.dumps(record))
    assert current.nitro_state_directory('g','s',policy)==folder
    record['model_identity']['reference_sha256']['atom_refs']='other refs'
    (folder/'result.json').write_text(json.dumps(record))
    with pytest.raises(ValueError):current.nitro_state_directory('g','s',policy)
