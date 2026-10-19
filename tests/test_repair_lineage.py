"""Canonical replacement must preserve the actual original source and SP identity."""
import hashlib

import pytest

from redox.qm import native_thermal_repair as repair
from redox.validation.provenance import repair_source


def test_promotion_source_resolves_after_canonical_geometry_is_replaced(tmp_path):
    canonical=tmp_path/'canonical';archive=tmp_path/'archive'
    canonical.mkdir();archive.mkdir()
    (archive/'opt.xyz').write_bytes(b'original geometry')
    expected=hashlib.sha256((archive/'opt.xyz').read_bytes()).hexdigest()
    (canonical/'opt.xyz').write_bytes(b'replacement geometry')
    snapshot=tmp_path/'producer.py';snapshot.write_bytes(repair.REPAIR_SOURCE_BYTES)
    original=dict(e_smd_eV=-123.5,g_thermal_eV=1.2,
                  minimum_repair_source=str(canonical/'opt.xyz'),
                  energy_protocol={'old':'producer'},energy_geom_sha1='old geometry')
    result=repair.promotion_metadata(original,canonical,archive,expected,snapshot)
    assert result['minimum_repair_source']==str(archive/'opt.xyz')
    assert hashlib.sha256((archive/'opt.xyz').read_bytes()).hexdigest()==result['minimum_repair_source_sha256']
    assert result['opt_out']==str(canonical/'opt.xyz')
    assert result['e_smd_eV']==original['e_smd_eV'] and result['g_thermal_eV']==original['g_thermal_eV']
    assert 'energy_protocol' not in result and 'energy_geom_sha1' not in result
    assert original['energy_geom_sha1']=='old geometry'  # private evidence is untouched


@pytest.mark.parametrize('corruption',['archive','producer'])
def test_promotion_refuses_unverified_archive_or_source(tmp_path,corruption):
    archive=tmp_path/'archive';archive.mkdir();p=archive/'opt.xyz';p.write_bytes(b'old')
    expected=hashlib.sha256(p.read_bytes()).hexdigest()
    snapshot=tmp_path/'producer.py';snapshot.write_bytes(repair.REPAIR_SOURCE_BYTES)
    (p if corruption=='archive' else snapshot).write_bytes(b'changed')
    with pytest.raises(RuntimeError,match='changed'):
        repair.promotion_metadata({},tmp_path/'canonical',archive,expected,snapshot)


def test_legacy_source_resolves_only_through_the_explicit_matching_archive(tmp_path):
    canonical=tmp_path/'canonical';archive=tmp_path/'archive';canonical.mkdir();archive.mkdir()
    (canonical/'opt.xyz').write_bytes(b'new geometry')
    p=archive/'opt.xyz';p.write_bytes(b'original geometry')
    raw=dict(minimum_repair_source='canonical/opt.xyz',
             minimum_repair_source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
             legacy_thermal_archive='archive/result.json')
    result=repair_source(raw,tmp_path)
    assert result['repair_source_status']=='resolved_archived_source'
    assert result['repair_source_geometry']==str(p)
    assert raw['minimum_repair_source']=='canonical/opt.xyz'
    p.write_bytes(b'different original')
    assert repair_source(raw,tmp_path)['repair_source_geometry'] is None


def test_source_resolution_never_falls_back_to_directory_names_or_missing_hashes(tmp_path):
    old=tmp_path/'superseded_canonical';old.mkdir();p=old/'opt.xyz';p.write_bytes(b'old')
    raw=dict(minimum_repair_source='missing.xyz',minimum_repair_source_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    assert repair_source(raw,tmp_path)['repair_source_status']=='unverified'
    raw['legacy_thermal_archive']='superseded_canonical/result.json';raw['minimum_repair_source_sha256']=None
    assert repair_source(raw,tmp_path)['repair_source_geometry'] is None
    assert repair_source({},tmp_path)['repair_source_status']=='not_declared'


def test_malformed_source_metadata_cannot_crash_status_collection(tmp_path):
    assert repair_source([],tmp_path)['repair_source_status']=='unverified'
    assert repair_source(None,tmp_path)['repair_source_geometry'] is None
