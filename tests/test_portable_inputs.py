import hashlib
import json

import pytest

from redox.validation.portable_inputs import PortableInputs


def inputs(tmp_path):
    historical=tmp_path/'original.py';historical.write_text('old producer')
    digest=hashlib.sha256(historical.read_bytes()).hexdigest()
    blob=tmp_path/'saved.py';blob.write_bytes(historical.read_bytes())
    manifest=tmp_path/'inputs.json'
    data=dict(schema=1,analysis_sha256='a'*64,bindings={str(historical):dict(file=blob.name,sha256=digest)})
    manifest.write_text(json.dumps(data))
    return historical,blob,manifest,data,digest


def test_historical_and_current_producers_remain_distinct(tmp_path):
    original,blob,path,data,old=inputs(tmp_path)
    reader=PortableInputs(path,{str(original):old},'a'*64)
    original.write_text('new producer');new=hashlib.sha256(original.read_bytes()).hexdigest()
    assert reader.verify(original,old)==blob
    with pytest.raises(ValueError):reader.verify(original,new)
    reader.allow_local([original])
    assert reader.verify(original,new)==original
    assert reader.verify(original,old)==blob


@pytest.mark.parametrize('change',['bytes','missing','traversal','symlink','analysis'])
def test_invalid_portable_binding_rejected(tmp_path,change):
    original,blob,path,data,digest=inputs(tmp_path)
    if change=='bytes':blob.write_text('corrupted')
    if change=='missing':data['bindings']={}
    if change=='traversal':data['bindings'][str(original)]['file']='../escape'
    if change=='symlink':blob.unlink();blob.symlink_to(original)
    if change=='analysis':data['analysis_sha256']='b'*64
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        reader=PortableInputs(path,{str(original):digest},'a'*64);reader.verify(original,digest)
