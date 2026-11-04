#!/usr/bin/env python3
"""Run the frozen polymer gradient diagnostic with an isolated screening variant."""
import argparse
from contextlib import nullcontext
import importlib.util
import json
from pathlib import Path
import runpy


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    args = parser.parse_args()
    from redox.core.io import atomic_write_json
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    case = args.plan.parent
    sources = case.parent / 'sources'
    plan = json.loads(args.plan.read_text())
    variant = plan['screening_variant']
    if variant not in ('installed', 'corrected'):
        raise ValueError('An explicit frozen screening variant is required')
    definition = importlib.util.spec_from_file_location('polymer_screening_helper', sources / 'helper.py')
    helper = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(helper)
    context = nullcontext() if variant == 'installed' else helper.corrected_screening(
        sources / 'kernel.cu', case / 'screening_implementation')
    atomic_write_json(case / 'screening_variant.json', dict(variant=variant,
        native_G_qualified=False, production_default_changed=False))
    with context:
        runpy.run_path(str(sources / 'scientific_worker.py'), run_name='__main__')
