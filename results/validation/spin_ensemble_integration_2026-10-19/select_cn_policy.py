"""Select only after the actual hook and its regression suite pass."""
import hashlib,json,os,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from redox.validation.spin_ensemble import select_policy
base=Path(__file__).resolve().parent
activation=json.loads((base/'activation.json').read_text());assert activation['status']=='activated_pending_full_checks'
log=ROOT/'results/reviews/2026-10-19/software_checks_spin_ensemble_activated.log';assert '422 passed' in log.read_text()
assert hashlib.sha256((ROOT/'src/redox/core/common.py').read_bytes()).hexdigest()==activation['candidate_sha256']
policy=json.loads((base/'prepared_policy.json').read_text())
result=select_policy(ROOT/'calcs/dft/search_pmdi_1864c2c74267/red2',policy)
result.update(tests_passed=422,test_log_sha256=hashlib.sha256(log.read_bytes()).hexdigest(),selected=time.time(),quantum_jobs_submitted=0)
atomic_write_json(base/'selection_receipt.json',result)
activation.update(status='activated_and_selected_readers_verified',software_tests_passed=422,selection_receipt=str(base/'selection_receipt.json'))
atomic_write_json(base/'activation.json',activation);print(result)
