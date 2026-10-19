"""Every thermally included spin component must pass integrity, not just its carrier."""
import importlib.util
from pathlib import Path

import pytest

from redox.core.common import ROOT
from redox.properties import integrity
from redox.validation import spin_ensemble


@pytest.mark.parametrize('alternative_energy,bond_changes,expected',[
    (-11.8,0,'ok'),(-9.,0,'unbound'),(-11.8,2,'bond_change'),(-11.8,None,'incomplete'),
])
def test_alternative_surface_controls_the_necessary_condition(tmp_path,monkeypatch,
                                                            alternative_energy,bond_changes,expected):
    ox=tmp_path/'m/neu/opt.xyz';red=tmp_path/'m/red1/opt.xyz';alt=tmp_path/'alternative/opt.xyz'
    for p in (ox,red,alt):p.parent.mkdir(parents=True)
    (red.parent/spin_ensemble.POLICY_NAME).write_text('{}')
    monkeypatch.setattr(spin_ensemble,'admitted_ensemble',lambda *args:dict(components=[
        dict(directory=str(red.parent),electronic_eV=-12.),
        dict(directory=str(alt.parent),electronic_eV=alternative_energy)]))
    monkeypatch.setattr(integrity,'connectivity_change',lambda a,b:bond_changes if b==alt else 0)
    rO,rR={'e_smd_eV':-10.},{'e_smd_eV':-12.}
    result=integrity.ensemble_integrity(rO,rR,ox,red,-1)
    assert result['ensemble_component_pairs']==2
    assert result['ensemble_integrity_qc']==expected
    monkeypatch.setattr(integrity,'DFT',tmp_path)
    monkeypatch.setattr(integrity,'_couples',lambda gid:[(('neu',0),('red1',-1))])
    monkeypatch.setattr(integrity,'_res',lambda gid,st:rO if st=='neu' else rR)
    monkeypatch.setattr(integrity,'kabsch_rmsd_heavy',lambda *args:0.)
    row,=integrity.assess('m')
    assert row['verdict']==('intact_bound' if expected=='ok' else expected)


def test_invalid_selected_ensemble_never_falls_back_to_one_geometry(tmp_path,monkeypatch):
    ox=tmp_path/'o/opt.xyz';red=tmp_path/'r/opt.xyz';red.parent.mkdir()
    (red.parent/spin_ensemble.POLICY_NAME).write_text('{}')
    monkeypatch.setattr(spin_ensemble,'admitted_ensemble',lambda *args:None)
    result=integrity.ensemble_integrity({'e_smd_eV':-10.},{'e_smd_eV':-12.},ox,red,-1)
    assert result['ensemble_integrity_qc']=='incomplete'


def test_independent_distance_graph_reconstructs_bond_breaking(tmp_path):
    spec=importlib.util.spec_from_file_location('independent_integrity_test',ROOT/'scripts/validation/audit/recompute_axes.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    short,long=tmp_path/'short.xyz',tmp_path/'long.xyz'
    short.write_text('2\n\nH 0 0 0\nH 0 0 .74\n')
    long.write_text('2\n\nH 0 0 0\nH 0 0 2.\n')
    assert audit.independent_bond_changes(short,short)==0
    assert audit.independent_bond_changes(short,long)==1
    assert integrity.connectivity_change(short,long)==1
