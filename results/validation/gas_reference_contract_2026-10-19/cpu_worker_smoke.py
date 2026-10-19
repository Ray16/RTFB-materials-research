"""Actual small CPU four-point calculation, saved-point audit and qualified resume."""
import hashlib,importlib.util,json,os,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

base=Path(__file__).resolve().parent;out=base/'cpu_worker_smoke';out.mkdir(exist_ok=False)
with allocation_guard():
 resources=host_capacity(1);cpus=cpu_allocation(1);os.sched_setaffinity(0,cpus)
assert all(os.environ[k]=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'))
source=ROOT/'scripts/validation/reorg_literature/validate_reorg_worker_generic.py'
spec=importlib.util.spec_from_file_location('actual_cpu_reference_worker',source);worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
worker.CALC=out/'calc';worker.XC='pbe';worker.BASIS='sto-3g'
row=dict(id='water_hole_control',smiles='O',lit_hole_eV=None)
started=time.time();status=worker.process(row,'hole','cpu');assert status=='ok',status
record=json.loads((worker.CALC/'water_hole_control.json').read_text());e=record['four_point_evidence']['energies_Ha']
independent=((e['ion_at_neutral']-e['ion_at_ion'])+(e['neutral_at_ion']-e['neutral_at_neutral']))*27.211386245988
assert abs(independent-record['our_lambda'])<=.00005001
assert record['solver_verification']['identity']['basis']=='sto-3g'
assert worker.process(row,'hole','cpu')=='skip'
report=dict(status='actual_CPU_reference_worker_and_resume_passed',system='water neutral/cation, PBE/STO-3G, density-fitted gas four-point',
 lambda_eV=record['our_lambda'],independent_unrounded_lambda_eV=independent,saved_four_points_verified=True,
 qualified_resume_skipped_recomputation=True,CPU_affinity=cpus,resources=resources,GPU_used=False,
 experimental_accuracy_validated=False,elapsed_seconds=time.time()-started,
 source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
 (Path(__file__),source,ROOT/'src/redox/qm/gas_reference.py',worker.CALC/'water_hole_control.json')})
atomic_write_json(out/'report.json',report);print(json.dumps(report,indent=2))
