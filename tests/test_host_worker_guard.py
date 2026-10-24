"""A host-memory guard preserves a successful worker and stops only itself."""
import json
from pathlib import Path
import subprocess
import sys


def command(tmp_path, limit, source):
    script = tmp_path/'worker.py'
    script.write_text(source)
    out = tmp_path/'guard'
    guard = Path(__file__).resolve().parents[1]/'scripts/fleet/bounded_host_python_worker.py'
    result = subprocess.run([sys.executable, str(guard), '--out', str(out),
        '--maximum-rss-gib', str(limit), '--minimum-free-gib', '.001', str(script)],
        capture_output=True, text=True, timeout=15)
    return result, json.loads((out/'status.json').read_text())


def test_completed_python_program_retains_same_pid_and_needs_review(tmp_path):
    marker = tmp_path/'pid.json'
    result, report = command(tmp_path, 1, f'import os,json\nfrom pathlib import Path\nPath({str(marker)!r}).write_text(json.dumps(os.getpid()))\n')
    assert result.returncode == 0
    assert report['pid'] == json.loads(marker.read_text())
    assert report['status'] == 'program_returned_review_required'
    assert report['peak_RSS_bytes'] > 0


def test_exceeded_bound_interrupts_worker_and_preserves_guard_record(tmp_path):
    result, report = command(tmp_path, .001, 'import time\ntime.sleep(10)\n')
    assert result.returncode != 0
    assert report['peak_RSS_bytes'] > report['maximum_RSS_GiB']*2**30
    assert report['signal_target'] == 'this owned process only'
    assert report['status'] == 'program_interrupted_or_failed_preserved'
    assert report['finished'] >= report['started']
