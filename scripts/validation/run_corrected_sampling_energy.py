#!/usr/bin/env python3
"""Run unchanged smooth-model gas/SMD hybrid equations under corrected screening.

The only storage adaptation is an explicit occupied-only determinant checkpoint;
all virtual coefficients remain available to the live SCF calculation. No old
electronic/nuclear value or incomplete polymer cycle is promoted to a result.
"""
import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time


def load(path, name):
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition); definition.loader.exec_module(module)
    return module


@contextmanager
def occupied_storage(api, driver_path, recorder):
    original = api.load_module
    def module(name, path):
        result = original(name, path)
        if Path(path) == driver_path:
            result.scf_checkpoints = recorder.scf_checkpoints
        return result
    api.load_module = module
    try:
        yield
    finally:
        api.load_module = original


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--plan', type=Path, required=True)
    a = p.parse_args(); case = a.plan.parent; sources = case.parent / 'sources'
    from redox.core.common import ROOT
    from redox.core.io import atomic_write_json
    from redox.workflows import polymer_campaign as api
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner(); plan = api.read_plan(case)
    if (plan['role'] != 'corrected_sampling_energy' or not plan['corrected_sampling_energy']
            or plan['checkpoint_format'] != 'redox.restricted_occupied_SCF.v1'
            or plan['historical_energy_or_thermal_values_transferred']
            or plan['native_G_computed'] or plan['primary_candidate_admission']):
        raise ValueError('Exact corrected electronic-only method-control scope required')
    def verify():
        if any(hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest
               for path, digest in plan['source_sha256'].items()):
            raise ValueError('Corrected sampling source snapshot changed')
    verify()
    record = dict(status='corrected_sampling_energy_running', started=time.time(),
        owner=dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
        native_G_qualified=False, occupied_ensemble_cycle_qualified=False)
    atomic_write_json(case / 'status.json', record)
    try:
        helper = load(sources / 'helper.py', 'sampling_corrected_screening')
        recorder = load(sources / 'occupied.py', 'sampling_occupied_checkpoint')
        driver = load(ROOT / 'scripts/pipeline/smooth_polymer_energy_bounded.py', 'sampling_original_bounded_energy')
        with helper.corrected_screening(sources / 'kernel.cu', case / 'screening_implementation'), \
                occupied_storage(api, ROOT / 'scripts/pipeline/smooth_polymer_energy.py', recorder):
            driver.run(case)
        verify()
        stage = json.loads((case / 'experimental_energy/status.json').read_text())
        if stage['status'] != 'experimental_hybrid_energy_finished_independent_review_pending':
            raise ValueError('Original gas/SMD energy execution did not finish normally')
        record['observed_peak_GPU_MiB'] = stage['observed_peak_GPU_MiB']
        record['status'] = 'corrected_sampling_energy_finished_independent_review_pending'
    except BaseException as error:
        record.update(status='corrected_sampling_energy_failed_preserved', error=repr(error)); raise
    finally:
        record['finished'] = time.time(); atomic_write_json(case / 'status.json', record)


if __name__ == '__main__':
    main()
