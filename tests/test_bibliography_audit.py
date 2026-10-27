import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('bibliography_audit',ROOT/'scripts/validation/audit_bibliography.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def test_nested_tex_and_unicode_metadata_match():
    row=m.entries(r'''@article{x, author={Mardirossian, Narbe and Others}, title={{$\omega$B97M-V}: {Gaussian} Tests}, doi={10.example/a},year={2016}}''')[0]
    record={'registered':dict(doi='10.EXAMPLE/a',title='ωB97M-V: Gaussian Tests',year_candidates=[2015,2016],first_author='Mardirossian')}
    assert m.compare(row,record)['status']=='metadata_matches'


@pytest.mark.parametrize('key,value',[('doi','wrong'),('title','Another paper'),('author','Other, Name'),('year','2000')])
def test_wrong_paper_identity_is_never_verified(key,value):
    row=dict(doi='10.x/y',title='Exact title',author='Smith, A',year='2020');row[key]=value
    record={'registered':dict(doi='10.x/y',title='Exact title',first_author='Smith',year_candidates=[2020])}
    assert m.compare(row,record)['status']=='manual_review_required'


def test_duplicate_keys_and_unclosed_fields_are_rejected():
    with pytest.raises(ValueError):m.entries('@misc{x,title={a}}\n@misc{x,title={b}}')
    with pytest.raises(ValueError):m.field('title={unclosed','title')


def test_tex_accent_commands_are_not_part_of_author_identity():
    assert m.normalized(r"Bari{\v s}i{\' c}")==m.normalized('Barišić')
    assert m.normalized(r'$\beta$')==m.normalized('β')


def test_relative_metadata_cache_replays_without_network(tmp_path,monkeypatch):
    import hashlib,json
    bib=tmp_path/'refs.bib';bib.write_text('@article{x,doi={10.x/y},title={Title},author={Smith, A},year={2020}}')
    cache=tmp_path/'old';cache.mkdir();record=dict(registered=dict(doi='10.x/y',title='Title',first_author='Smith',year_candidates=[2020]))
    metadata=cache/'metadata.json';metadata.write_text(json.dumps(record));summary=cache/'summary.json'
    summary.write_text(json.dumps(dict(checks=[dict(status='metadata_matches',citation={'doi':'10.x/y'},metadata_file='metadata.json',metadata_sha256=hashlib.sha256(metadata.read_bytes()).hexdigest())])))
    monkeypatch.setattr(m,'get_metadata',lambda doi:pytest.fail('Must reuse exact cache without network'))
    result=m.audit(bib,tmp_path/'new',summary)
    assert result['checks'][0]['status']=='metadata_matches'


@pytest.mark.parametrize('key',['../outside','a/b'])
def test_citation_keys_cannot_escape_metadata_directory(key):
    with pytest.raises(ValueError,match='citation key'):
        m.entries('@misc{'+key+',title={Title}}')
