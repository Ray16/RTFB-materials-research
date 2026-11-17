import copy

import pytest

from redox.workflows.ionic_stationary_cycles import electronic_cycle


@pytest.fixture
def cycle():
    # A small balanced chemical bookkeeping example: ethane + H2 - 2 methane.
    definition=dict(coefficients=dict(both=1,empty=1,first=-1,second=-1))
    cases={key:dict(smiles=smiles,charge=0) for key,smiles in [('both','[H][H]'),('empty','CC'),('first','C'),('second','C')]}
    records={key:dict(E_Ha=value,G_eV=999.) for key,value in [('both',-5.),('empty',-1.),('first',-2.),('second',-3.)]}
    return definition,records,cases


def test_balanced_cycle_reports_electronic_difference_but_never_promotes_G(cycle):
    result=electronic_cycle(*cycle)
    assert result['delta_delta_E_Ha']==-1.
    assert result['delta_delta_G_kJ_mol'] is None
    assert not result['isolated_motif_binding_energy'] and not result['electronic_coupling_computed']


def test_missing_component_remains_pending(cycle):
    definition,records,cases=cycle;records.pop('empty');result=electronic_cycle(definition,records,cases)
    assert result['missing']==['empty'] and result['delta_delta_E_Ha'] is None and result['delta_delta_G_kJ_mol'] is None


@pytest.mark.parametrize('change',['atoms','charge','coefficient','nonfinite'])
def test_bad_stoichiometry_or_nonfinite_cycle_cannot_pass(cycle,change):
    definition,records,cases=copy.deepcopy(cycle)
    if change=='atoms':cases['empty']['smiles']='CCC'
    if change=='charge':cases['both']['charge']=1
    if change=='coefficient':definition['coefficients']['both']=.5
    if change=='nonfinite':records['both']['E_Ha']=float('nan')
    with pytest.raises(ValueError):electronic_cycle(definition,records,cases)


def test_cycle_implementation_receipt_must_bind_actual_corrected_kernel(tmp_path):
    import hashlib
    from redox.workflows.ionic_stationary_cycles import screening_signature
    kernel=tmp_path/'kernel.cu';compiled=tmp_path/'screening_kernel.cu';kernel.write_bytes(b'kernel');compiled.write_bytes(b'patched kernel')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt=dict(status='run_local_screening_context_closed_qualification_separate',upstream_fix='https://github.com/pyscf/gpu4pyscf/pull/925',
        kernel_sha256=sha(compiled),source_sha256=sha(kernel),installed_gen_grid_sha256='installed',installed_library_files_modified=False)
    assert screening_signature(receipt,tmp_path,compiled)['kernel_sha256']==sha(compiled)
    kernel.write_bytes(b'changed')
    with pytest.raises(ValueError):screening_signature(receipt,tmp_path,compiled)
