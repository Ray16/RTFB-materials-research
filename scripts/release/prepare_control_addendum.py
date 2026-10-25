#!/usr/bin/env python
"""Package checked raw model controls and a path-independent CPU replay.

No quantum calculation, record repair, model substitution or source-state write
is performed. A public archival license/DOI remains a separate author decision.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def function_copy(source, names, target, imports):
    text = source.read_text(); tree = ast.parse(text)
    functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    if not set(names) <= set(functions):
        raise ValueError('Expected independent audit functions are missing')
    target.write_text('"""Frozen independent replay functions; provenance is in the manifest."""\n'
        + imports + '\n\n' + '\n\n'.join(ast.get_source_segment(text, functions[name]) for name in names)+'\n')


def prepare(out, numpy_lock):
    panels = [('CHO', 'smooth_quantum_2026-10-25', {1:'formaldehyde', 2:'methyl'}),
              ('NCl', 'smooth_quantum_NCl_2026-10-25', {1:'ammonia', 2:'chloromethyl'})]
    frozen = ROOT / 'manuscript/evidence/2026-10-25_reference_diagnostic'
    frozen_manifest = json.loads((frozen / 'manifest.json').read_text())
    api.verify_sources({str(ROOT / p):v for p,v in frozen_manifest['snapshot_sha256'].items()})
    reviews = {}
    for family, name, labels in panels:
        path = ROOT / 'results/validation' / name / 'independent_review.json'
        review = json.loads(path.read_text())
        if review['status'] != 'experimental_smooth_quantum_controls_independently_checked':
            raise ValueError('All molecular controls require completed independent source reviews')
        api.verify_sources(review['source_sha256']); reviews[family] = review
    fc_root = ROOT / 'results/validation/fc_functional_sensitivity_2026-10-25'
    fc_review = json.loads((fc_root / 'independent_review.json').read_text())
    api.verify_sources(fc_review['source_sha256'])
    out.mkdir(parents=True, exist_ok=False)
    example = out / 'example'; example.mkdir()
    data = example / 'data'; data.mkdir()
    original = {}
    def copy(source, target):
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        expected = digest(source)
        if digest(target) != expected: raise ValueError('Copied raw bytes changed')
        original[str(source)] = expected
        return str(target.relative_to(example))
    copy(ROOT / 'examples/control_replay/reproduce.py', example / 'reproduce.py')
    copy(numpy_lock, example / 'requirements-cpu.lock.txt')
    copies = [(ROOT / 'manuscript/check_claims.py', ['audit_fc_functional_controls', 'audit_surface_controls'],
               'claims.py', 'import math'),
              (ROOT / 'scripts/validation/review_smooth_quantum.py', ['projected_modes'],
               'projected_modes.py', 'import numpy as np'),
              (ROOT / 'scripts/validation/audit/native_partition.py', ['nuclear_equations'],
               'nuclear_equations.py', 'import math\nimport numpy as np')]
    derived = {}
    for source, functions, name, imports in copies:
        function_copy(source, functions, example / name, imports)
        original[str(source)] = digest(source)
        derived[name] = dict(original_source=str(source), source_sha256=digest(source), functions=functions)
    evidence_names = ['Fc_functional_review.json', 'benchmark.json', 'CDS_contact.json',
        'CDS_smooth_component.json', 'CDS_smooth_qualification.json', 'CDS_smooth_full_gradients.json',
        'CDS_smooth_quantum_CHO.json', 'CDS_smooth_quantum_NCl.json']
    evidence_names += [f'Fc_functional_{method}_{state}.json'
        for method in ('wb97mv', 'b3lyp', 'pbe0') for state in ('neu', 'ox')]
    evidence = {name:copy(frozen / name, data / 'evidence' / name) for name in evidence_names}
    quantum = []
    for family, name, labels in panels:
        root = ROOT / 'results/validation' / name
        for mult, species in labels.items():
            source = json.loads((root / 'controls' / str(mult) / 'status.json').read_text())
            if source['status'] != 'two_model_quantum_force_Hessian_and_nuclear_controls_passed':
                raise ValueError('Source quantum control did not finish normally')
            api.verify_sources(source['output_sha256'])
            # Keep runtime/ownership records as immutable descriptive evidence.
            for role, path in [('source',root/'controls'/str(mult)/'status.json'),
                ('runtime',root/'runtime'/str(mult)/'status.json'),
                ('guard_owner',root/'guards'/str(mult)/'owner.json'),
                ('guard_state',root/'guards'/str(mult)/'state.json')]:
                copy(path, data / 'provenance' / f'{species}_{role}.json')
            report = next(r for r in reviews[family]['controls'] if r['mult'] == mult)
            for arm in ('legacy', 'smooth'):
                raw = Path(source['raw_root']) / arm
                target = data / 'quantum' / (species+'__'+arm)
                fields = {'geometry':'opt.xyz', 'energy':'quantum_result.json', 'nuclear':'nuclear.json',
                          'hessian':'hessian.npy', 'gradient':'gradient.npy'}
                files = {role:copy(raw / filename, target / filename) for role,filename in fields.items()}
                evaluations = sorted((raw / 'SCF_evaluations').glob('*/evaluation.json'))
                if len(evaluations) != 49:
                    raise ValueError('The raw small-control displacement set is incomplete')
                points = [copy(p, target / 'evaluations' / (p.parent.name+'.json')) for p in evaluations]
                quantum.append(dict(species=species, arm=arm, files=files, evaluations=points,
                    reference_metrics=report['arms'][arm]))
    fc = []
    for method in ('wb97mv', 'b3lyp', 'pbe0'):
        for state in ('neu', 'ox'):
            case = fc_root / (method+'__'+state)
            geometry = copy(case / 'input.xyz', data / 'Fc' / (method+'__'+state+'.xyz'))
            for role, path in [('source',case/'electronic/status.json'),
                ('runtime',case/'tensor_runtime/status.json'),('guard_owner',case/'residency_guard/owner.json'),
                ('guard_state',case/'residency_guard/state.json')]:
                copy(path, data / 'provenance' / f'Fc_{method}_{state}_{role}.json')
            fc.append(dict(method=method, state=state, geometry=geometry,
                energy_evidence=f'Fc_functional_{method}_{state}.json'))
    (example / 'README.md').write_text('''# Cached Fc and solvent-surface control replay

This small CPU example contains original bytes for six Fc electronic records
and eight small quantum control states (four molecules, two surface models).
Each quantum state includes geometry, active gas/SMD energies, a Cartesian DFT
Hessian/gradient, nuclear data and all 49 baseline/displaced SCF records.

Standard-library electronic/component replay, from any directory:

`python3 -I -S /path/to/example/reproduce.py --electronic-only`

For full cached derivative/mode/nuclear replay, create a Python 3.11 Linux x86_64
venv, install the hash-locked NumPy wheel, then run:

`python -m pip install --require-hashes -r /path/to/example/requirements-cpu.lock.txt`

`python -I /path/to/example/reproduce.py`

The replay checks checksums, contained paths, exact state/geometry/model joins,
Fc electronic gaps/counterfactuals, component curvatures, twelve force components
and two Hessian columns at two displacements per control, and independently
reconstructed modes and nuclear sums at three rotor cutoffs. The comparison
controls impose rotational symmetry number one; this is not the primary-state
symmetry procedure. No original project/calculation path is opened.

This replays stored numerical evidence. It does not run SCF/DFT, load CUDA or
xTB, repeat the physical calculations, validate predictive solvent accuracy,
admit a candidate, establish polymer G or demonstrate ensemble convergence.
Code is MIT-licensed. Final public data-deposition licensing/DOI is pending.
''')
    model = json.loads((frozen / 'CDS_smooth_qualification.json').read_text())['method']
    manifest = dict(schema=1, primary_candidate_admission=False, polymer_free_energy_qualification=False,
        experimental_model=model, quantum_states=quantum, electronic_states=fc, evidence=evidence,
        derived_source_functions=derived, source_evidence_epoch=str(frozen.relative_to(ROOT)),
        source_evidence_manifest_sha256=digest(frozen / 'manifest.json'),
        scope='Cached control arithmetic and raw derivative/mode/nuclear reconstruction; no new quantum or primary admission',
        files_sha256={str(p.relative_to(example)):digest(p) for p in example.rglob('*') if p.is_file()})
    (data / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    # Preserve the exact already-public source revision, plus the new replay's
    # own source bytes separately; neither reads the user's unrelated index.
    commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with (out / 'code_snapshot.tar').open('xb') as stream:
        subprocess.run(['git','archive','--format=tar',commit,'src','scripts','config','environments',
            'pyproject.toml','README.md','LICENSE'],cwd=ROOT,stdout=stream,check=True)
    copy_source = out / 'prepare_source.py'; shutil.copyfile(Path(__file__).resolve(),copy_source)
    environment = out / 'environments'; environment.mkdir()
    for path in (ROOT / 'environments').glob('*'):
        if path.is_file() and (path.name.endswith('.txt') or path.suffix == '.json'):
            shutil.copyfile(path, environment / path.name)
    (environment / 'cutensor22').mkdir()
    shutil.copyfile(ROOT/'environments/cutensor22/runtime.json',environment/'cutensor22/runtime.json')
    api.verify_sources(original)
    provenance = dict(status='raw_control_addendum_prepared_relocation_pending', prepared=time.time(),
        source_code_commit=commit, source_sha256=original, quantum_control_state_records=8,
        Fc_electronic_state_records=6, raw_SCF_evaluations=8*49,baseline_evaluations=8,displaced_evaluations=8*48,
        original_records_modified=False, quantum_jobs_submitted=0,
        primary_candidate_admission=False, public_archive_deposited=False,
        scope='Complement to the v5 and first-site bundles; overlapping identities are not new unique molecules')
    (out / 'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (out / 'README.md').write_text('''# Raw numerical-control addendum

Run the self-contained CPU example under `example/`. The outer package also
contains the exact public code revision, recorded quantum environment artifacts
and raw-byte source provenance. The cached CPU replay does not install or
qualify the full GPU environment. The example's NumPy requirement is separately
hash-locked for Python 3.11 Linux x86_64. CUDA/model weights are not redistributed.

This package complements, rather than replaces, the earlier v5 and first-site
packages. It supplies no new polymer G or converged ensemble. Public archival
metadata, data license and DOI remain pending author confirmation.
''')
    package = dict(schema=1, status='control_addendum_files_frozen_relocation_pending',
        source_code_commit=commit, files_sha256={str(p.relative_to(out)):digest(p)
            for p in out.rglob('*') if p.is_file()})
    (out / 'package_manifest.json').write_text(json.dumps(package,indent=2)+'\n')
    print(json.dumps(dict(status=package['status'], package=str(out), files=len(package['files_sha256']),
        quantum_control_states=8, electronic_states=6, raw_SCF_evaluations=392,baseline_evaluations=8,displaced_evaluations=384)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--numpy-lock', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard(): host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1))
    prepare(args.out.resolve(), args.numpy_lock.resolve())


if __name__ == '__main__': main()
