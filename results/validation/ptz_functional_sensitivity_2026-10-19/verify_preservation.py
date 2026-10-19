"""Compare primary numerical artifacts with the preceding published checkpoint."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json

base = Path(__file__).resolve().parent
checkpoint = '075cc61b1479c7256babb0fe97cd621ce05fa7e9'
tables = 'results/reviews/2026-10-19/uncertainty_scope_tables'
paths = subprocess.check_output(
    ['git', 'ls-tree', '-r', '--name-only', checkpoint, '--', tables],
    cwd=ROOT, text=True).splitlines()
paths = [p for p in paths if p.endswith(('.csv', '.xlsx'))]
paths += ['results/reviews/2026-10-18/verified_candidates/LATEST_verified_candidates.xlsx',
          'results/reviews/2026-10-18/verified_candidates/latest.json',
          'config/model_validation_limits.json', 'results/validated/current',
          'results/candidates/latest.xlsx', 'results/candidates/latest.csv',
          'results/candidates/comparisons.csv', 'results/candidates/current.json']
checks = []
for name in paths:
    expected = subprocess.check_output(['git', 'show', checkpoint + ':' + name], cwd=ROOT)
    path = ROOT / name
    actual = os.readlink(path).encode() if path.is_symlink() else path.read_bytes()
    checks.append(dict(path=name, sha256=hashlib.sha256(actual).hexdigest(),
                       unchanged=actual == expected))
relaxation = base.parent / 'ptz_active_relaxation_2026-10-19'
prior_plan = json.loads((relaxation / 'plan.json').read_text())
for task in prior_plan['tasks']:
    for name, expected in task['input_sha256'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == expected
    selected = json.loads((Path(task['source']) / 'selected_SP.json').read_text())
    geometry = Path(selected['geom'])
    result = geometry.parent / 'result.json'
    for path, expected in ((geometry, hashlib.sha256((Path(task['source']) / 'geometry.xyz').read_bytes()).hexdigest()),
                           (result, selected['adopted_result_sha256'])):
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        checks.append(dict(path=str(path.relative_to(ROOT)), sha256=actual,
                           unchanged=actual == expected, baseline='previously published frozen parent inputs'))
assert checks and all(c['unchanged'] for c in checks), [c for c in checks if not c['unchanged']]
baseline = ROOT / 'results/reviews/2026-10-19/git_cached_before_checkpoint.patch'
assert subprocess.check_output(['git', 'diff', '--cached', '--binary'], cwd=ROOT) == baseline.read_bytes()
report = dict(status='canonical_parent_candidate_and_property_tables_unchanged',
              baseline_checkpoint=checkpoint, checks=checks, failed_checks=0,
              unrelated_staged_patch_unchanged=True, checked=time.time())
atomic_write_json(base / 'preservation_review.json', report)
print(f"Unchanged: {len(checks)} primary artifacts; unrelated staged patch preserved")
