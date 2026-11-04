"""Maintained QC helpers need explicit snapshots, not archival path fallback."""
from pathlib import Path
import shutil

import pytest

from redox.core.hashing import sha256_file
from redox.validation.molecular_native import MolecularSources


def test_maintained_helper_snapshot_is_unambiguous_and_archive_reads_stay_strict(tmp_path):
    prior=tmp_path/'current';archive=tmp_path/'archive';case=tmp_path/'case'
    archive.mkdir();(prior/'src').mkdir(parents=True);(case/'review_sources').mkdir(parents=True)
    current=prior/'src/new_density_helper.py';current.write_text('maintained QC source')
    copied=case/'review_sources/density_helper.py';shutil.copyfile(current,copied)
    reader=MolecularSources(dict(prior_root=str(prior),archive=str(archive),case=str(case)))
    assert reader.verify(str(copied),sha256_file(copied))==copied
    with pytest.raises(FileNotFoundError):reader.verify(str(current),sha256_file(current))
    current.write_text('later maintained revision')
    assert reader.verify(str(copied),sha256_file(copied))==copied
