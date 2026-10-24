"""Collect unreachable SCF objects before creating the next GPU mean field.

Live scanners and densities remain reachable and are unchanged. This composes
the qualified workspace adapter; full spin controls qualify the new lifecycle.
"""
from contextlib import contextmanager
import gc
from pathlib import Path
import time


@contextmanager
def collect_before_build(report):
    from redox.qm import dft
    previous=dft._build_smd_mf
    def build(*args,**kwargs):
        collected=gc.collect()
        dft._free_gpu_pool()
        report['collections']+=1
        report['unreachable_objects_collected']+=collected
        return previous(*args,**kwargs)
    dft._build_smd_mf=build
    try:yield
    finally:dft._build_smd_mf=previous


@contextmanager
def compact_smd_execution(directory,budget_GiB):
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api
    parent=ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py'
    module=api.load_module('GC_workspace_parent',parent)
    report=dict(status='GC_workspace_active_pending_controls',started=time.time(),collections=0,
        unreachable_objects_collected=0,scientific_model_changed=False,
        source_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(),parent)})
    with module.compact_smd_execution(directory,budget_GiB) as evidence,collect_before_build(report):
        try:
            yield evidence
            report['status']='GC_workspace_finished_requires_control_and_physical_checks'
        except BaseException as error:
            report.update(status='GC_workspace_failed_preserved',error=repr(error));raise
        finally:
            report['finished']=time.time();atomic_write_json(Path(directory)/'GC_execution.json',report)
