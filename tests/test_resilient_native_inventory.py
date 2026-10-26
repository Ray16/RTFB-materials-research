import importlib.util
from pathlib import Path
import subprocess
from types import SimpleNamespace

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('inventory_adapter',ROOT/'scripts/pipeline/coordinate_native_inventory_resilient.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


def test_inventory_reports_verified_slots_and_unreachable_hosts_separately():
    def fake(cmd,**kwargs):
        if cmd[-1]=='lambda1':return SimpleNamespace(stdout='lambda1 2\nlambda1 3\n')
        raise subprocess.TimeoutExpired(cmd,15)
    a=mod.query_host('lambda1',run=fake);b=mod.query_host('lambda2',run=fake)
    assert a['status']=='verified' and a['slots']==[('lambda1',2),('lambda1',3)]
    assert b['status']=='inventory_unverified' and not b['slots']


def test_inventory_rejects_misattribution_and_duplicate_slots():
    for raw in ['lambda2 0\n','lambda1 0\nlambda1 0\n','invalid\n']:
        result=mod.query_host('lambda1',run=lambda *a,**k:SimpleNamespace(stdout=raw))
        assert result['status']=='inventory_unverified' and not result['slots']
