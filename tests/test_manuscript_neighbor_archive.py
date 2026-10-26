from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('neighbor_archive_claims_test',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def evidence():
    folder=ROOT/'manuscript/evidence/2026-10-26_native_replay'
    return {name:json.loads((folder/name).read_text()) for name in
        ['neighbor_archive.json','neighbor_archive_relocation.json','neighbor_archive_manifest.json','neighbor_cycle.json']}


def test_complete_raw_cycle_and_native_replay_support_the_manuscript(evidence):
    result=module.audit_neighbor_archive(evidence.__getitem__)
    assert result['Hessians_reconstructed']==13 and result['raw_states']==14
    assert result['independent_nuclear_reconstruction'] and not result['original_quantum_QC_repeated']


@pytest.mark.parametrize('change',['archive','counts','original_read','NumPy','NaN','state_G','scope','cycle','missing_Hessian'])
def test_neighbor_archive_claims_reject_mismatched_or_promoted_evidence(evidence,change):
    evidence=deepcopy(evidence);proof=evidence['neighbor_archive_relocation.json']
    if change=='archive':proof['archive_sha256']='another-archive'
    elif change=='counts':proof['archive_regular_files']+=1
    elif change=='original_read':proof['native_isolation']['attempted_original_reads']=['original.xyz']
    elif change=='NumPy':proof['native_isolation']['numpy_version']='unqualified'
    elif change=='NaN':proof['native_reconstruction']['states'][0]['nuclear_errors_eV']['50']=float('nan')
    elif change=='state_G':next(r for r in evidence['neighbor_archive_manifest.json']['states'] if r['id']=='empty')['checked_G_eV']+=1
    elif change=='scope':evidence['neighbor_archive_manifest.json']['conformer_convergence_established']=True
    elif change=='cycle':proof['arithmetic']['by_cutoff']['50']['deltaDeltaG_kJ_mol']+=1
    else:proof['native_reconstruction']['states'].pop()
    with pytest.raises(ValueError):module.audit_neighbor_archive(evidence.__getitem__)
