"""Collection may reclaim dead cycles but must preserve a live SCF object."""
import importlib.util
from pathlib import Path
import weakref

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('gc_workspace',ROOT/'scripts/pipeline/smd_gc_workspace.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def test_collects_only_unreachable_cycles_and_restores_builder(monkeypatch):
    from redox.qm import dft
    class Object:pass
    dead=Object();dead.cycle=dead;dead_ref=weakref.ref(dead);del dead
    live=Object();live.cycle=live;live_ref=weakref.ref(live)
    calls=[]
    def builder(*args,**kwargs):calls.append((args,kwargs));return live
    monkeypatch.setattr(dft,'_build_smd_mf',builder);monkeypatch.setattr(dft,'_free_gpu_pool',lambda:None)
    proof=dict(collections=0,unreachable_objects_collected=0)
    with MODULE.collect_before_build(proof):
        assert dft._build_smd_mf('molecule',xc='unchanged') is live
        assert dead_ref() is None and live_ref() is live
    assert dft._build_smd_mf is builder and calls==[(('molecule',),dict(xc='unchanged'))]
    assert proof['collections']==1
