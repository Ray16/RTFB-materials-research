import importlib.util
import io
from pathlib import Path
import tarfile

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('validate_manuscript_bundle',ROOT/'scripts/release/validate_manuscript_bundle.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def archive(entries):
    out=io.BytesIO()
    with tarfile.open(fileobj=out,mode='w') as stream:
        for name,kind in entries:
            item=tarfile.TarInfo(name)
            if kind=='symlink':item.type=tarfile.SYMTYPE;item.linkname='../outside';stream.addfile(item)
            else:
                content=b'data';item.size=len(content);stream.addfile(item,io.BytesIO(content))
    out.seek(0);return out


def test_safe_extraction_rejects_escape_symlink_and_duplicate(tmp_path):
    for index,entries in enumerate([[('../escape','file')],[('root/link','symlink')],[('root/a','file'),('root/a','file')]]):
        out=tmp_path/str(index);out.mkdir()
        with tarfile.open(fileobj=archive(entries),mode='r:') as source,pytest.raises((ValueError,FileExistsError)):
            module.safe_extract(source,out,expected_root='root')
    assert not (tmp_path/'escape').exists()


def test_safe_extraction_preserves_declared_regular_file(tmp_path):
    with tarfile.open(fileobj=archive([('root/a','file')]),mode='r:') as source:
        files,roots=module.safe_extract(source,tmp_path,expected_root='root')
    assert files=={'root/a'} and roots=={'root'} and (tmp_path/'root/a').read_bytes()==b'data'
