"""Run the five CPU replays against captured, successful source-tree baselines."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def main():
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1',
                       OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                       PYTHONPATH=str(ROOT / 'src'))
    failures = []
    for manifest in ('runs.json', 'scientific_runs.json'):
        for record in json.loads((ROOT / 'baseline' / manifest).read_text()):
            if record['returncode'] != 0:
                continue
            command = [sys.executable, *record['command'][1:]]
            expected = record.get('stdout_sha256', record.get('sha256'))
            try:
                result = subprocess.run(command, cwd=ROOT, env=environment,
                                        capture_output=True, timeout=180)
                actual = hashlib.sha256(result.stdout).hexdigest()
                passed = result.returncode == 0 and actual == expected
                report = dict(command=command, returncode=result.returncode,
                              expected=expected, actual=actual, passed=passed)
                if not passed:
                    report['stderr'] = result.stderr.decode(errors='replace')
            except subprocess.TimeoutExpired:
                passed = False
                report = dict(command=command, passed=False, error='timeout after 180s')
            print(json.dumps(report), flush=True)
            if not passed:
                failures.append(command)
    return int(bool(failures))


if __name__ == '__main__':
    raise SystemExit(main())
