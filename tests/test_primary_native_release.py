"""A native resource pilot cannot silently admit a different scientific class."""
import json
from pathlib import Path

import pytest

from redox.core.hashing import sha256_file
from redox.workflows import screening_native_panel as panel

QUANTUM=('native_thermal.py','projected_thermal.py','qrrho.py','thermal.py','dft.py','cds_finite.py','hessian_memory.py','xtb_kernel.py','vibrations.py')


def setup(tmp_path,monkeypatch):
    archive=tmp_path/'archive';archive.mkdir();witness=tmp_path/'physical_receipt';witness.write_text('unit-test receipt')
    specs={};rows=[];reviews=[]
    for key,mult,nao in [('rks_pilot',1,394),('uks_pilot',2,394),('rks_child',1,175),('uks_child',2,221)]:
        root=tmp_path/key;case=root/'case';sources=root/'sources';case.mkdir(parents=True);sources.mkdir()
        for name in ('worker.py','helper.py','checkpoints.py','hashing.py','kernel.cu'):(sources/name).write_text('same frozen test source')
        plan=dict(role='primary_molecular_Hessian_control',native_nao=nao,mult=mult,charge=1 if mult==2 else 2,
                  native_memory_cap_GiB=6,source_sha256={str(archive/'src/redox/qm'/n):'same-bounded-source' for n in QUANTUM})
        (case/'plan.json').write_text(json.dumps(plan));(case/'result.json').write_text(json.dumps(dict(G_by_cutoff_eV={c:-100. for c in ('25','50','100')})))
        ns=root/'namespace.json';ns.write_text('{}');specs[str(ns)]=dict(case=str(case),sources=str(sources),archive=str(archive))
        rows.append(dict(key=key,mult=mult,nao=nao,namespace=str(ns),namespace_sha256=sha256_file(ns)))
        if key.endswith('pilot'):
            report=dict(status='primary_molecular_native_control_independently_reconstructed',key=key,namespace=str(ns),
                charge=plan['charge'],mult=mult,native_component_requalified=True,material_admission=False,
                conformer_ensemble_qualified=False,actual_SCF_checkpoint_frames_checked=9,peak_GPU_MiB=23000,
                maximum_column_error=1e-5,minimum_frequency_cm=10.,source_sha256={str(witness):sha256_file(witness)},
                G_by_cutoff_eV={c:-100. for c in ('25','50','100')})
            p=tmp_path/(key+'.json');p.write_text(json.dumps(report));reviews.append(p)
    path=tmp_path/'manifest.json';path.write_text(json.dumps(dict(states=rows,target_native_states=len(rows),pilot_keys=['rks_pilot','uks_pilot'])))
    monkeypatch.setattr(panel,'verify',lambda p:specs[str(p)])
    return path,reviews


def test_exact_completed_RKS_UKS_pilots_release_only_their_shared_class(tmp_path,monkeypatch):
    manifest,reviews=setup(tmp_path,monkeypatch)
    assert len(panel.qualified_rows(manifest,reviews))==4


@pytest.mark.parametrize('bad',['missing_pilot','wrong_G','wrong_spin','nonfinite','worker','dimension','thermal','pending'])
def test_missing_or_inconsistent_pilots_and_changed_implementations_stay_held(tmp_path,monkeypatch,bad):
    manifest,reviews=setup(tmp_path,monkeypatch)
    if bad=='missing_pilot':reviews.pop()
    if bad in ('wrong_G','wrong_spin','nonfinite','pending'):
        p=reviews[0];r=json.loads(p.read_text())
        if bad=='wrong_G':r['G_by_cutoff_eV']['50']+=.01
        if bad=='wrong_spin':r['mult']=3
        if bad=='nonfinite':r['G_by_cutoff_eV']['25']=float('nan')
        if bad=='pending':r['native_component_requalified']=False
        p.write_text(json.dumps(r))
    if bad=='worker':(tmp_path/'rks_child/sources/worker.py').write_text('changed')
    if bad in ('dimension','thermal'):
        p=tmp_path/'rks_child/case/plan.json';r=json.loads(p.read_text())
        if bad=='dimension':r['native_nao']=395
        if bad=='thermal':r['source_sha256'][str(tmp_path/'archive/src/redox/qm/qrrho.py')]='different source'
        p.write_text(json.dumps(r))
    with pytest.raises(ValueError):panel.qualified_rows(manifest,reviews)
