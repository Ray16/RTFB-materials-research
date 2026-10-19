"""Experimental calibration requires current computed endpoints and spin decisions."""
import pytest

from redox.core import common
from redox.properties import potentials
from redox.validation import reference_admission as gate


@pytest.fixture
def admitted(monkeypatch):
    records={key:dict(e_smd_eV=value,g_thermal_eV=0,thermal_qc='ok',
                      thermal_definition='native_solution_dft_qrrho_spin')
             for key,value in [(('m','ox'),5.5),(('m','neu'),1.),(('ferrocene','ox'),4.)]}
    monkeypatch.setattr(common,'read_result',lambda gid,state:records.get((gid,state)))
    monkeypatch.setattr(potentials,'fc_live',lambda:4.)
    monkeypatch.setattr(gate,'spin_comparison_running',lambda gid:False)
    monkeypatch.setattr(gate,'state_spin_evidence',lambda *a:dict(spin_primary_admissible=True))
    return records


def test_matching_current_event_is_allowed(admitted):
    assert gate.event_admission('m','ox->neu',.5)[0]


@pytest.mark.parametrize('value',[.501,float('nan'),float('inf'),True])
def test_stale_or_nonfinite_serialization_is_rejected(admitted,value):
    assert not gate.event_admission('m','ox->neu',value)[0]


def test_missing_native_endpoint_is_rejected(admitted):
    admitted['m','neu']['g_thermal_eV']=None
    assert not gate.event_admission('m','ox->neu',.5)[0]


def test_unmatched_thermal_reference_is_rejected(admitted):
    admitted['ferrocene','ox']['thermal_definition']='legacy_gfn2_gas_rrho'
    assert not gate.event_admission('m','ox->neu',.5)[0]


def test_pending_spin_endpoint_is_rejected(admitted,monkeypatch):
    monkeypatch.setattr(gate,'state_spin_evidence',lambda *a:
                        dict(spin_primary_admissible=False,spin_confidence='matching native spin thermals pending'))
    assert not gate.event_admission('m','ox->neu',.5)[0]


def test_active_spin_writer_blocks_calibration(admitted,monkeypatch):
    monkeypatch.setattr(gate,'spin_comparison_running',lambda gid:True)
    assert not gate.event_admission('m','ox->neu',.5)[0]


def test_grounded_measurement_cannot_enter_sigma_without_computed_admission(monkeypatch):
    from redox.screening import scorecard
    monkeypatch.setattr(gate,'event_admission',lambda *a:(False,'unqualified selected spin'))
    assert scorecard.experimental_e_points([
        dict(id='tempo_parent',event='ox->rad',E_vs_Fc_V='.449')])==[]
