"""Full-contact queue admission preserves the pilot's physical/resource scope."""
from copy import deepcopy
import json

import pytest

from redox.core.hashing import sha256_file
from redox.workflows import counterion_thermochemistry as api


def fixture(tmp_path,monkeypatch):
    case=tmp_path/'case';(case/'state').mkdir(parents=True);xyz=case/'state/opt.xyz';xyz.write_text('checked exact coordinates')
    ns=tmp_path/'namespace.json';ns.write_text('{}')
    pilot=dict(electronic_protocol={'xc':'primary'},native_protocol={'Hessian':'primary'},
               thermal_protocol={'T':298.15},smiles='exact contact',peak_host_RSS_bytes=2**30)
    plan=dict(role='primary_counterion_thermochemistry',pilot_review='checked pilot',charge=1,mult=1,atoms=42,
              dimensions={'def2-svp':412,'def2-tzvpd':1003},native_G_qualified=False,
              maximum_peak_GPU_MiB=31744,maximum_host_RSS_GiB=16,geometry_sha256=sha256_file(xyz),
              **{k:deepcopy(pilot[k]) for k in ('electronic_protocol','native_protocol','thermal_protocol','smiles')})
    manifest=tmp_path/'manifest.json';manifest.write_text(json.dumps(dict(states=[dict(key='N0_below_plus',namespace=str(ns),namespace_sha256=sha256_file(ns))])))
    monkeypatch.setattr(api,'verify',lambda p:dict(case=str(case)))
    monkeypatch.setattr(api,'checked_pilot',lambda p:(pilot,{}))
    return case,manifest,plan


def test_queue_admits_only_prepared_full_class_and_supplies_explicit_launcher(tmp_path,monkeypatch):
    case,manifest,plan=fixture(tmp_path,monkeypatch);(case/'plan.json').write_text(json.dumps(plan))
    assert api.prepared_rows(manifest)[0]['launcher_module']=='redox.workflows.counterion_thermochemistry'


@pytest.mark.parametrize('field,value',[('minimum_GPU_memory_MiB',15360),('maximum_peak_GPU_MiB',15360),
    ('mult',2),('charge',0),('dimensions',{'def2-svp':397}),('pilot_review',None),
    ('thermal_protocol',{'T':300}),('native_G_qualified',True),('maximum_host_RSS_GiB',8)])
def test_small_card_radical_context_and_unqualified_resource_changes_are_refused(tmp_path,monkeypatch,field,value):
    case,manifest,plan=fixture(tmp_path,monkeypatch);plan[field]=value;(case/'plan.json').write_text(json.dumps(plan))
    with pytest.raises(ValueError):api.prepared_rows(manifest)
