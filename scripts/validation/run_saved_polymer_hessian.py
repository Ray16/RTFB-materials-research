#!/usr/bin/env python3
"""Run the unchanged physical checks using a closed complete Hessian checkpoint."""
import argparse
import importlib.util
import json
from pathlib import Path


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);a=p.parse_args()
    plan=json.loads(a.plan.read_text());case=a.plan.parent;sources=case.parent/'sources'
    if not plan.get('Hessian_validation_continuation') or not plan.get('saved_full_Hessian'):
        raise ValueError('Explicit source-bound saved-Hessian validation plan required')
    from redox.qm import dft
    helper=load(sources/'saved_matrix_helper.py','exact_saved_Hessian_builder')
    scientific=load(sources/'scientific_worker.py','unchanged_full_Hessian_physical_checks')
    with helper.builder(dft,plan,case):scientific.main()


if __name__=='__main__':main()
