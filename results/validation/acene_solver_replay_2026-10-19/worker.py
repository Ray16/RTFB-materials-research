"""One fresh gas hole-relaxation anchor; leave every legacy byte intact."""
import argparse,hashlib,importlib.util,json,os,socket,time
from pathlib import Path
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner

base=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--task',type=int,required=True);args=parser.parse_args()
plan=json.loads((base/'plan.json').read_text());task=plan['tasks'][args.task]
folder=base/'tasks'/task['row']['id'];folder.mkdir(parents=True,exist_ok=False)
report=dict(status='starting',task=task,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),experimental_accuracy_validated=False,production_promotion=False)
def save():atomic_write_json(folder/'status.json',report)
def verify():
 for p,h in {**plan['producer_sha256'],**task['input_sha256']}.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
save()
try:
 verify_gate_owner();verify()
 assert len(os.sched_getaffinity(0))==2
 assert all(os.environ[k]=='2' for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'])
 spec=importlib.util.spec_from_file_location('current_acene_worker',plan['worker_source']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 module.CALC=base/'calc';module.BASIS='6-31g(d,p)'
 report.update(status='computing',physical_gpu=os.environ['CUDA_VISIBLE_DEVICES']);save()
 status=module.process(task['row'],'gpu');assert status=='ok',status
 record=json.loads((module.CALC/(task['row']['id']+'.json')).read_text());assert record['status']=='ok' and record['solver_verification']['contract']['full_grid_response']
 assert len(record['four_point_evidence']['energies_Ha'])==4
 verify();report.update(status='fresh_solver_qualified_record_written',result=str(module.CALC/(task['row']['id']+'.json')),result_sha256=hashlib.sha256((module.CALC/(task['row']['id']+'.json')).read_bytes()).hexdigest(),our_hole_eV=record['our_hole'])
except BaseException as exc:
 report.update(status='failed_preserving_outputs',error=repr(exc));raise
finally:
 report['finished']=time.time();save()
