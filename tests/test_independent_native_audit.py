import hashlib
import importlib.util
import json

import pytest

from redox.core.common import ROOT


def audit_module():
    spec = importlib.util.spec_from_file_location('native_audit_test', ROOT / 'scripts/validation/audit/recompute_axes.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture_record(tmp_path, audit):
    xyz = '2\n\nH 0 0 0\nH 0 0 .7\n'
    (tmp_path/'opt.xyz').write_text(xyz)
    definition = 'native_solution_dft_qrrho_spin'
    (tmp_path/'thermal_policy.json').write_text(json.dumps(dict(thermal_definition=definition,
                                                               native_protocol=audit.NATIVE_THERMAL_CONTRACT)))
    ident = dict(geometry_sha256=hashlib.sha256(xyz.encode()).hexdigest(), charge=0, mult=1,
                 protocol=audit.NATIVE_THERMAL_CONTRACT)
    key = hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
    path = tmp_path/'native_thermal'/f'{key}.json';path.parent.mkdir()
    proof=[];controls=[]
    for index,status in enumerate(('analytic_full_response_verified','uks_full_response_verified',
                                  'fine_response_implementation_verified')):
        control=tmp_path/f'control{index}.json';control.write_text(json.dumps(dict(status=status)))
        proof.append(dict(path=str(control),status=status,sha256=hashlib.sha256(control.read_bytes()).hexdigest()))
        controls.append(dict(relative_path=str(control),source_sha256=proof[-1]['sha256'],source_status=status,
                             admission_status=status,software=dict(pyscf='2.14.0',gpu4pyscf='1.8.1'),
                             relative_full_matrix_error=.0001,fine_response_implementation_verified=True))
        if index==2:proof[-1]['relative_full_matrix_error']=.0001
    certificate=dict(schema=1,native_protocol=audit.NATIVE_THERMAL_CONTRACT,
                     software=dict(pyscf='2.14.0',gpu4pyscf='1.8.1'),controls=controls)
    certificate['certificate_content_sha256']=hashlib.sha256(json.dumps(certificate,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    audit.NATIVE_CERTIFICATE=tmp_path/'certificate.json';audit.NATIVE_CERTIFICATE.write_text(json.dumps(certificate))
    record=dict(identity=ident,status='derivative_checked_solution_thermal',converged=True,
                stationary_minimum=True,grid_response=True,frequencies_cm=[4000.],gradient_max=.0001,
                maximum_asymmetry=1e-10,derivative_checks=[dict(atom_index=i,axis=axis,
                step_Bohr=step,relative_error=.0001) for i,axis in ((0,1),(1,2)) for step in (.0005,.001)],
                smd_spin_diagnostics=dict(spin_contam=0.,homo_eV=-5.,s_squared=0.),
                software=dict(pyscf='2.14.0',gpu4pyscf='1.8.1'),implementation_controls=proof,
                g_nuclear_qrrho_eV=1.,g_spin_partition_eV=0.,g_thermal_total_eV=1.,nuclear_by_cutoff_eV={'50':1.})
    path.write_text(json.dumps(record))
    return path,record,dict(charge=0,mult=1,g_thermal_eV=99.,n_imag=0)


def test_audit_admits_selected_native_data_without_legacy_fallback(tmp_path):
    audit=audit_module();path,record,raw=fixture_record(tmp_path,audit)
    assert audit._surface_thermal_inputs(tmp_path,raw)==(1.,0,'native_solution_dft_qrrho_spin')
    path.unlink()
    assert audit._surface_thermal_inputs(tmp_path,raw)[0] is None


def test_independent_primary_requires_corrected_partition(tmp_path,monkeypatch):
    audit=audit_module();path,record,raw=fixture_record(tmp_path,audit)
    assert audit._surface_thermal_inputs(tmp_path,raw)[0]==1.
    assert audit.thermal_inputs(tmp_path,raw)[0] is None
    monkeypatch.setattr(audit,'native_partition_inputs',lambda *args:dict(g_thermal_total_eV=.75))
    assert audit.thermal_inputs(tmp_path,raw)==(.75,0,'native_solution_dft_qrrho_spin')
    record['maximum_asymmetry']=.01;path.write_text(json.dumps(record))
    assert audit.thermal_inputs(tmp_path,raw)[0] is None


@pytest.mark.parametrize('change',[dict(g_thermal_total_eV=2.),dict(frequencies_cm=[-20.]),
                                 dict(derivative_checks=[]),dict(maximum_asymmetry=float('nan'))])
def test_audit_independently_rejects_bad_native_records(tmp_path,change):
    audit=audit_module();path,record,raw=fixture_record(tmp_path,audit)
    record.update(change);path.write_text(json.dumps(record))
    assert audit._surface_thermal_inputs(tmp_path,raw)[0] is None
