import importlib.util
import io
from pathlib import Path
import sys
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
definition = importlib.util.spec_from_file_location('native_cache_isolation',
    ROOT/'scripts/release/validate_native_cache_isolation.py')
validator = importlib.util.module_from_spec(definition)
sys.path.insert(0, str(ROOT/'scripts/release'))
try:
    definition.loader.exec_module(validator)
finally:
    sys.path.pop(0)


def stream(names, kind=tarfile.REGTYPE):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode='w') as archive:
        for name in names:
            item = tarfile.TarInfo(name); item.type = kind
            item.size = 4 if kind == tarfile.REGTYPE else 0
            archive.addfile(item, io.BytesIO(b'data') if item.size else None)
    data.seek(0)
    return tarfile.open(fileobj=data, mode='r|')


def test_expected_regular_members_extract_to_the_declared_root(tmp_path):
    with stream(['package/manifest.json', 'package/objects/hash']) as archive:
        assert validator.safe_extract(archive, tmp_path, 'package') == 2
    assert (tmp_path/'package/objects/hash').read_bytes() == b'data'


@pytest.mark.parametrize('name', ['/outside', '../outside', 'package/../../outside', 'other/record'])
def test_unexpected_roots_and_parent_paths_are_rejected(tmp_path, name):
    with stream([name]) as archive, pytest.raises(ValueError):
        validator.safe_extract(archive, tmp_path, 'package')


@pytest.mark.parametrize('kind', [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.CHRTYPE])
def test_links_and_devices_cannot_reach_original_files(tmp_path, kind):
    with stream(['package/record'], kind) as archive, pytest.raises(ValueError):
        validator.safe_extract(archive, tmp_path, 'package')


def test_duplicate_archive_files_are_not_overwritten(tmp_path):
    with stream(['package/record', 'package/record']) as archive, pytest.raises(ValueError, match='Duplicate'):
        validator.safe_extract(archive, tmp_path, 'package')


def test_previous_archive_is_preserved_before_any_packaging(tmp_path):
    archive = tmp_path/'existing.tar.zst'
    archive.write_bytes(b'previous evidence')
    with pytest.raises(FileExistsError):
        validator.compress(tmp_path/'absent_package', archive)
    assert archive.read_bytes() == b'previous evidence'
