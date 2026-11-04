import copy
from pathlib import Path

import pytest

from redox.validation.workspace_index import relocate_index


def fixture(tmp_path):
    root=tmp_path/'project';new=root/'calcs'/'case';new.mkdir(parents=True)
    (new/'value.json').write_text('{}')
    old=tmp_path/'old_case';old.symlink_to(new,target_is_directory=True)
    mapping={'moved':[{'old':str(old),'new':str(new)}]}
    return root,old,new,mapping


def test_exact_pointer_changes_preserve_status_and_deferred_locations(tmp_path):
    root,old,new,mapping=fixture(tmp_path)
    original={'raw_root':str(old),'record':str(old/'value.json'),'scope':'old_case is historical',
              'pending':None,'QC':{'accepted':False},'live':str(tmp_path/'live_campaign')}
    saved=copy.deepcopy(original)
    result,changes=relocate_index(original,mapping,root)
    assert original==saved and len(changes)==2
    assert result==dict(saved,raw_root=str(new),record=str(new/'value.json'))
    assert relocate_index(result,mapping,root)==(result,[])


@pytest.mark.parametrize('damage',['missing','outside','foreign_old','traversal','overlap'])
def test_invalid_relocation_cannot_change_the_index(tmp_path,damage):
    root,old,new,mapping=fixture(tmp_path)
    source={'record':str(old/'value.json')};saved=copy.deepcopy(source)
    if damage=='missing':source['record']=str(old/'absent.json');saved=copy.deepcopy(source)
    if damage=='outside':mapping['moved'][0]['new']=str(tmp_path)
    if damage=='foreign_old':old.unlink();old.mkdir()
    if damage=='traversal':source['record']=str(old/'..'/'escape');saved=copy.deepcopy(source)
    if damage=='overlap':mapping['moved'].append(dict(mapping['moved'][0]))
    with pytest.raises(ValueError):relocate_index(source,mapping,root)
    assert source==saved


def test_prefix_lookalikes_are_not_relocated(tmp_path):
    root,old,new,mapping=fixture(tmp_path)
    source={'path':str(old)+'_another/file'}
    assert relocate_index(source,mapping,root)==(source,[])


def test_input_lookup_works_after_the_compatibility_link_is_removed(tmp_path):
    from redox.validation.workspace_index import relocated_input_path
    root,old,new,mapping=fixture(tmp_path)
    old.unlink()
    assert relocated_input_path(old/'value.json',mapping,root)==new/'value.json'
    assert (new/'value.json').read_text()=='{}'
    assert relocated_input_path(tmp_path/'unmapped.json',mapping,root)==tmp_path/'unmapped.json'
