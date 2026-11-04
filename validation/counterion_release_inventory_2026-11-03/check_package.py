"""Verify staged objects with original scientific paths blocked."""
import json
import os
from pathlib import Path
import subprocess
import sys

from redox.core.hashing import sha256_file
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

with allocation_guard():
    host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT.parent / 'redox_release_candidates/2026-11-03_declared_contact_sources'
OUT = Path(__file__).with_name('package_verification.json')
if OUT.exists():
    raise FileExistsError('Preserve previous verification')
blocked = [ROOT, ROOT.parent/'archive', ROOT.parent/'redox_counterion_controls',
           ROOT.parent/'redox_candidate_continuations', ROOT.parent/'redox_kernel_controls']
code = '''import json,pathlib,runpy,sys
package=pathlib.Path(sys.argv[1]);blocked=[pathlib.Path(x).resolve() for x in sys.argv[2:]]
def deny(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes)):
  path=pathlib.Path(args[0]).resolve()
  if any(path==p or path.is_relative_to(p) for p in blocked):raise PermissionError('Original scientific paths blocked')
sys.addaudithook(deny)
namespace=runpy.run_path(str(package/'verify.py'),run_name='isolated_verifier')
print(json.dumps(namespace['verify'](package)))
'''
assert sha256_file(PACKAGE/'verify.py') == sha256_file(ROOT/'scripts/release/verify_source_objects.py')
run = subprocess.run([sys.executable, '-I', '-S', '-c', code, str(PACKAGE),
                      *map(str, blocked)], cwd='/', text=True, capture_output=True, check=True)
result = json.loads(run.stdout)
result.update(package=str(PACKAGE), manifest_sha256=sha256_file(PACKAGE/'manifest.json'),
              verifier_sha256=sha256_file(PACKAGE/'verify.py'),
              original_scientific_paths_blocked=True,
              producer_sha256=sha256_file(__file__))
OUT.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
