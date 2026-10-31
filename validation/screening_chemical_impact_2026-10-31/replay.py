#!/usr/bin/env python3
"""Reproduce the cached 71-state/144-column checks with standard Python only."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path


def equivalent(actual,expected):
    if isinstance(expected,dict):
        return isinstance(actual,dict) and set(actual)==set(expected) and all(equivalent(actual[k],v) for k,v in expected.items())
    if isinstance(expected,list):
        return isinstance(actual,list) and len(actual)==len(expected) and all(equivalent(a,b) for a,b in zip(actual,expected))
    if type(expected) is float:
        return type(actual) in (int,float) and math.isfinite(actual) and abs(actual-expected)<=1e-14+1e-12*abs(expected)
    return actual==expected


def replay(root):
    manifest=json.loads((root/'manifest.json').read_text())
    for name,wanted in manifest['files_sha256'].items():
        path=(root/name).resolve()
        if not path.is_relative_to(root.resolve()) or hashlib.sha256(path.read_bytes()).hexdigest()!=wanted:
            raise ValueError('Changed or escaped replay artifact: '+name)
    source=root/'sensitivity_replay_source.py'
    spec=importlib.util.spec_from_file_location('frozen_screening_math',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    energy=json.loads((root/'energy_inputs.json').read_text());molecule=json.loads((root/'molecular_inputs.json').read_text())
    actual_energy=module.reconstruct(energy['states'],energy['endpoints'],energy['hartree_eV'])
    actual_molecule=module.molecular_derivatives(molecule['cases'])
    if not equivalent(actual_energy,json.loads((root/'energy_replay.json').read_text())) or not equivalent(actual_molecule,json.loads((root/'molecular_replay.json').read_text())):
        raise ValueError('Saved electronic/derivative evidence does not reconstruct')
    return dict(status='stdlib_screening_evidence_replay_passed',electronic_states=actual_energy['states'],
        endpoints=actual_energy['endpoints'],derivative_checks=actual_molecule['checks'],
        maximum_corrected_installed_electronic_shift_V=actual_energy['maximum_electronic_potential_shifts_V']['corrected_minus_installed_V'],
        quantum_calculations_run=False,native_G_qualified=False)


if __name__=='__main__':
    print(json.dumps(replay(Path(__file__).resolve().parent)))
