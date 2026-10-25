#!/usr/bin/env python
"""Bind the isolated fixed-orbital derivative probe to a qualified runtime."""
import argparse
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run'))
    for name in ('source','controls','runtime-controls','out'):
        parser.add_argument('--'+name,type=Path,required=name=='out')
    args=parser.parse_args()
    original=ROOT/'scripts/validation/qualify_occupied_df_probe.py'
    module=api.load_module('fixed_orbital_tensor_runtime_probe',original)
    if args.action=='run':
        module.run(args.out.resolve());return
    if not args.source or not args.controls or not args.runtime_controls:
        parser.error('prepare requires source and both control roots')
    checker_path=ROOT/'scripts/validation/review_cutensor_controls.py'
    checker=api.load_module('fixed_orbital_tensor_runtime_controls',checker_path)
    evidence=checker.checked_controls(args.runtime_controls.resolve())
    module.prepare(args.source.resolve(),args.controls.resolve(),32,args.out.resolve())
    plan=api.read_plan(args.out.resolve())
    comparison=args.runtime_controls.resolve()/'comparison.json'
    profile=Path(json.loads(comparison.read_text())['runtime_profile'])
    plan.update(runtime_profile=str(profile),runtime_controls=str(args.runtime_controls.resolve()))
    plan['producer_sha256'].update(evidence)
    for p in (Path(__file__).resolve(),original,checker_path,ROOT/'scripts/fleet/tensor_runtime_worker.py'):
        plan['producer_sha256'][str(p)]=api.digest(p)
    plan['input_sha256'][str(profile)]=api.digest(profile)
    atomic_write_json(args.out.resolve()/'plan.json',plan)


if __name__=='__main__':main()
