"""Activate only the staged analysis hook after both live source guards are past."""
import argparse,ast,copy,hashlib,json,os,shlex,subprocess,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
b=Path(__file__).resolve().parent;jobs_root=ROOT/'results/validation/ndi_methyl_radical_minimum_2026-10-19'
a=argparse.ArgumentParser();a.add_argument('--apply',action='store_true');args=a.parse_args()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
activation=json.loads((b/'activation.json').read_text());source=Path(activation['source']);candidate=b/'common_candidate.py'
assert sha(candidate)==activation['candidate_sha256']
old=ast.parse((b/'common_before.py').read_text());new=ast.parse(candidate.read_text());filtered=[]
for node in copy.deepcopy(new.body):
 if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)):
  if node.name=='apply_thermal_policy':continue
  if node.name=='_apply_single_surface_thermal_policy':node.name='apply_thermal_policy'
 filtered.append(node)
assert ast.dump(ast.Module(body=filtered,type_ignores=[]))==ast.dump(old),'Only the analysis-reader wrapper may change'
plan=json.loads((jobs_root/'source_plus.json').read_text());workers=json.loads((jobs_root/'launch.json').read_text())['workers']
assert plan['producer_sha256'][str(source)]==activation['expected_source_sha256']
remote=r'''import json,sys,psutil
from pathlib import Path
base=Path(sys.argv[1]);out=[]
for sign in ['minus','plus']:
 matches=[]
 for p in psutil.process_iter(['pid','cmdline']):
  argv=p.info['cmdline'] or []
  if str(base/'branch.py') in argv and '--branch' in argv and argv[argv.index('--branch')+1]==sign:
   matches.append(dict(pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv))
 out.append(dict(branch=sign,matching_processes=matches))
print(json.dumps(out))'''
py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
r=subprocess.run(['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=10','lambda13',shlex.join([py,'-c',remote,str(jobs_root)])],capture_output=True,text=True,timeout=30)
assert r.returncode==0,r.stderr
observed=json.loads(r.stdout);proof=[]
for entry in observed:
 sign=entry['branch'];case=jobs_root/sign/'qualification/candidate/status.json';passed=False
 if case.is_file():
  report=json.loads(case.read_text());execution_plan=json.loads((case.parent.parent/'plan.json').read_text())
  assert execution_plan['producer_sha256']==plan['producer_sha256']
  assert report['task']['source']==str(jobs_root/sign) and report['task']['id']==plan['id']
  # This trusted qualifier writes its first report only after all frozen-source
  # and task-input checks. Retain the exact report seen at activation time.
  assert sha(jobs_root/'qualify.py')==plan['producer_sha256'][str(jobs_root/'qualify.py')]
  assert report['status'] in ['running','checking_worst_contact','isolated_target_verified','failed_preserved']
  passed=True;entry['guard_report']=report
 entry['source_check_passed']=passed
 entry['safe_to_change_analysis_source']=passed or not entry['matching_processes']
 proof.append(entry)
ready=all(p['safe_to_change_analysis_source'] for p in proof)
observation=dict(status='source_guards_clear' if ready else 'waiting_for_source_guard',workers=proof,checked=time.time(),activation_applied=False)
atomic_write_json(b/'source_guard_observation.json',observation)
print(observation['status'])
if not args.apply or not ready:raise SystemExit(0)
if sha(source)==activation['candidate_sha256']:
 print('Reader hook already active');raise SystemExit(0)
assert sha(source)==activation['expected_source_sha256'],'Shared source changed; preserve and inspect'
for path,h in plan['producer_sha256'].items():assert sha(path)==h,path
assert sha(b/'common_before.py')==activation['expected_source_sha256']
tmp=source.with_name('.common.ensemble.tmp');tmp.write_bytes(candidate.read_bytes());os.replace(tmp,source)
activation.update(status='activated_pending_full_checks',activated=time.time(),source_guard_proof=proof,reader_wrapper_only_AST_verified=True)
atomic_write_json(b/'activation.json',activation);print('Activated shared-reader hook after verified source checks')
