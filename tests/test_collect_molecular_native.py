import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('tested_molecular_collection',ROOT/'scripts/pipeline/collect_molecular_native.py')
collector=importlib.util.module_from_spec(definition);definition.loader.exec_module(collector)


@pytest.fixture
def case(tmp_path,monkeypatch):
    archive=tmp_path/'archive';archive.mkdir();captured=tmp_path/'captured';captured.write_text('original source')
    namespaces=[];rows=[]
    for state,q in [('ox2',2),('ox1',1),('neu',0)]:
        namespace=tmp_path/(state+'_namespace.json')
        namespace.write_text(json.dumps(dict(prior_root='/old',archive=str(archive),case=str(tmp_path/state),state=state)))
        values={str(c):float(q) for c in (25,50,100)}
        review=tmp_path/(state+'_review.json');review.write_text(json.dumps(dict(status='exact_molecular_native_G_independently_checked',
            namespace=str(namespace),namespace_sha256=collector.sha256_file(namespace),state=state,charge=q,
            candidate_alias='exact',product_inchikey='same',G_by_cutoff_eV=values,
            source_sha256={str(captured):collector.sha256_file(captured)})))
        namespaces.append(str(namespace));rows.append(dict(namespace=str(namespace),state=state,
            status='exact_molecular_native_G_independently_checked',review=str(review),review_sha256=collector.sha256_file(review),G_by_cutoff_eV=values))
    prepared=tmp_path/'prepared.json';prepared.write_text(json.dumps(dict(namespaces=namespaces)))
    source=tmp_path/'collection.json';source.write_text(json.dumps(dict(entries=[*rows,
        dict(namespace='older_failed_neutral',state='neu',status='failed_native_attempt_preserved')])))
    monkeypatch.setattr(collector,'verify',lambda p:json.loads(p.read_text()))
    return prepared,source,rows,captured,tmp_path/'out.json'


def test_explicit_recovered_namespace_join_preserves_old_failure_history(case):
    prepared,source,rows,captured,out=case
    result=collector.collect(prepared,[source],out)
    assert result['checked']==3 and not result['older_failures_relabelled']
    assert any(r['status']=='failed_native_attempt_preserved' for r in result['original_outcome_history'])
    assert [r['namespace'] for r in result['entries']]==json.loads(prepared.read_text())['namespaces']


@pytest.mark.parametrize('change',['missing','duplicate','review_bytes','changed_source','different_graph'])
def test_ambiguous_incomplete_or_changed_join_cannot_be_published(case,change):
    prepared,source,rows,captured,out=case
    if change=='missing':source.write_text(json.dumps(dict(entries=rows[:-1])))
    if change=='duplicate':source.write_text(json.dumps(dict(entries=rows+[rows[0]])))
    if change=='review_bytes':Path(rows[0]['review']).write_text('{}')
    if change=='changed_source':captured.write_text('changed')
    if change=='different_graph':
        p=Path(rows[0]['review']);d=json.loads(p.read_text());d['product_inchikey']='another';p.write_text(json.dumps(d))
        rows[0]['review_sha256']=collector.sha256_file(p);source.write_text(json.dumps(dict(entries=rows)))
    with pytest.raises(ValueError):collector.collect(prepared,[source],out)
