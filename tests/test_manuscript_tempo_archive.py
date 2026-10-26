import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('tempo_archive_claim_audit',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def data():
    root=ROOT/'manuscript/evidence/2026-10-26_tempo_raw_replay'
    names=['tempo_archive.json','tempo_archive_relocation.json','tempo_archive_manifest.json','tempo_native_potential.json']
    records={n:json.loads((root/n).read_text()) for n in names}
    digest=lambda name:hashlib.sha256((root/name).read_bytes()).hexdigest()
    return records,digest


def test_archive_claim_is_reconstructed_without_publication_claim():
    records,digest=data();r=module.audit_tempo_archive(lambda n:copy.deepcopy(records[n]),digest)
    assert r['raw_states']==4 and r['Hessians_reconstructed']==4
    assert not r['public_archive_deposited'] and not r['quantum_calculations_repeated']


@pytest.mark.parametrize('mutation',['hessians','blocked_reads','accuracy','potential'])
def test_archive_claim_rejects_overstatement(mutation):
    records,digest=data();r=records['tempo_archive_relocation.json']
    if mutation=='hessians':r['native_result']['native_reconstruction']['Hessians_reconstructed']=3
    if mutation=='blocked_reads':r['stdlib_isolation']['original_project_open_attempts']=['forbidden']
    if mutation=='accuracy':r['native_result']['experimental_accuracy_inferred']=True
    if mutation=='potential':r['native_result']['oxidation_V_vs_Fc_by_cutoff']['50']+=.1
    with pytest.raises(ValueError):module.audit_tempo_archive(lambda n:copy.deepcopy(records[n]),digest)
