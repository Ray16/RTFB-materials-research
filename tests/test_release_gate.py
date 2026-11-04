"""Incomplete or mismatched verification receipts cannot finalize a package."""
import copy

import pytest

from redox.validation.release_gate import verified_archive_entry


def inputs(tmp_path):
    archive=tmp_path/'raw.tar.zst';archive.write_bytes(b'fixture archive bytes')
    manifest={'objects':{'a'*64:123},'support_files_sha256':{'verify.py':'b'*64},'inventory_sha256':'c'*64}
    inventory={'inventory_sha256':'c'*64,'archive_objects':1,'research_evidence_bytes':123}
    receipt={'status':'all_streamed_archive_members_verified','manifest_sha256':'d'*64,'objects':1,
             'members_verified':3,'archive':str(archive),'archive_sha256':'e'*64,
             'archive_bytes':archive.stat().st_size,'numerical_replay_performed':False}
    return receipt,manifest,inventory


def test_matching_byte_verification_receipt_can_be_attached(tmp_path):
    receipt,manifest,inventory=inputs(tmp_path)
    saved=copy.deepcopy((receipt,manifest,inventory))
    entry=verified_archive_entry(receipt,manifest,'d'*64,inventory,tmp_path)
    assert entry['file']=='raw.tar.zst' and entry['objects']==1
    assert (receipt,manifest,inventory)==saved


@pytest.mark.parametrize('damage',['pending','manifest','inventory','object_count','members','bytes','digest','replay','file','outside'])
def test_incomplete_or_wrong_receipt_is_refused(tmp_path,damage):
    receipt,manifest,inventory=inputs(tmp_path)
    if damage=='pending':receipt['status']='source_checked_archive_created_verification_pending'
    if damage=='manifest':receipt['manifest_sha256']='f'*64
    if damage=='inventory':inventory['inventory_sha256']='f'*64
    if damage=='object_count':receipt['objects']=2
    if damage=='members':receipt['members_verified']=2
    if damage=='bytes':inventory['research_evidence_bytes']=124
    if damage=='digest':receipt['archive_sha256']='not a digest'
    if damage=='replay':receipt['numerical_replay_performed']=True
    if damage=='file':(tmp_path/'raw.tar.zst').write_bytes(b'changed')
    if damage=='outside':receipt['archive']=str(tmp_path.parent/'external.tar.zst')
    with pytest.raises(ValueError):verified_archive_entry(receipt,manifest,'d'*64,inventory,tmp_path)
