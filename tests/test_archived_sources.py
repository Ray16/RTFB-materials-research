import hashlib
from pathlib import Path

import pytest

from redox.validation.archived_sources import ArchivedSources


def make(tmp_path):
    old, archive = tmp_path/'old', tmp_path/'archive'
    old.mkdir()
    archive.mkdir()
    return old, archive, ArchivedSources(old, archive)


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()


def test_explicit_relocation_does_not_trust_current_checkout(tmp_path):
    old, archive, reader = make(tmp_path)
    (old/'value').write_text('changed checkout')
    (archive/'value').write_text('frozen')
    assert reader.verify(old/'value', digest('frozen')) == archive/'value'
    assert (old/'value').read_text() == 'changed checkout'


def test_missing_archive_never_falls_back_to_old_path(tmp_path):
    old, _, reader = make(tmp_path)
    (old/'value').write_text('frozen')
    with pytest.raises(FileNotFoundError):
        reader.verify(old/'value', digest('frozen'))


def test_cache_rejects_rewritten_same_size_file(tmp_path):
    old, archive, reader = make(tmp_path)
    (archive/'value').write_text('first')
    reader.verify(old/'value', digest('first'))
    (archive/'value').write_text('other')
    with pytest.raises(ValueError, match='differs'):
        reader.recheck()


def test_external_paths_are_verified_without_prefix_substitution(tmp_path):
    old, _, reader = make(tmp_path)
    external = tmp_path/'old_extra'
    external.write_text('external')
    assert reader.verify(external, digest('external')) == external
    assert reader.path(external) != reader.archive/'_extra'


@pytest.mark.parametrize('path', ['relative/file', '/tmp/../input'])
def test_traversal_or_relative_paths_are_not_provenance(tmp_path, path):
    _, _, reader = make(tmp_path)
    with pytest.raises(ValueError):
        reader.path(path)


def test_digest_mismatch_is_not_a_relocation_exception(tmp_path):
    old, archive, reader = make(tmp_path)
    (archive/'value').write_text('altered')
    with pytest.raises(ValueError):
        reader.verify(old/'value', digest('frozen'))


def test_relative_bindings_require_explicit_original_directory(tmp_path):
    old, archive, _ = make(tmp_path)
    (archive/'value').write_text('frozen')
    reader = ArchivedSources(old, archive, relative_to_prior=True)
    assert reader.verify('value', digest('frozen')) == archive/'value'
    with pytest.raises(ValueError):
        reader.path('../outside')
