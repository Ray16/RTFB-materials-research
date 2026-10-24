from pathlib import Path
import json,time,fcntl,subprocess,hashlib,shutil
R=Path.cwd();root=R/'results/polymer_models/campaigns/2026-10-23_crest_refinement';audit=R/'results/validation/gpu_exclusivity_2026-10-23'
old=root/'baseline__design_pmdi_ethyl__occ_2_3_exclusive_resume';new=root/'baseline__design_pmdi_ethyl__occ_2_3_exclusive_second_resume';base=R/'results/validation/pmdi_occupancy_2026-10-23'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
 tmp=p.with_name(p.name+'.second_resume.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
events=[json.loads(s) for s in (audit/'lambda13/events.jsonl').read_text().splitlines() if s.strip()]
e=[e for e in events if e['pid']==2359770 and e['process_start']=='92391705'];assert len(e)==1
state=json.loads((old/'energy/workers/lambda13_gpu0_queue.json').read_text());assert state['status']=='stopped_preserving_outputs'
remote=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','lambda13',"python3 -c 'import json,pathlib,subprocess; p=pathlib.Path(\"/proc/2359770/stat\"); print(json.dumps(dict(stopped_worker_stat=p.read_text() if p.exists() else None,residency=subprocess.check_output([\"nvidia-smi\",\"--query-compute-apps=gpu_uuid,pid,used_gpu_memory\",\"--format=csv,noheader\"],text=True))))'"],text=True,capture_output=True,check=True)
observed=json.loads(remote.stdout);assert observed['stopped_worker_stat'] is None
proof=dict(checked=time.time(),event=e[0],subsequent_observation=observed,unrelated_processes_signaled=False,source_worker=str(old/'energy/workers/lambda13_gpu0_queue.json'),source_worker_sha256=sha(old/'energy/workers/lambda13_gpu0_queue.json'),scientific_failure_inferred=False,new_continuation=str(new),warm_start_source='Original complete evaluation 23; interrupted exclusive attempt had no new accepted result')
p=audit/'second_residency_interruption.json';write(p,proof)
plan=json.loads((new/'plan.json').read_text());plan['input_sha256'][str(p)]=sha(p);plan['prior_exclusive_attempt']=dict(campaign=str(old),stop_proof=str(p),no_accepted_E_or_G=True)
for k in ('producer_sha256','input_sha256'):
 for f,h in plan[k].items():assert sha(f)==h,(k,f)
write(new/'plan.json',plan)
readme=(old/'README.md').read_text();(old/'README.md').write_text(readme+'\nThe subsequent actual start was interrupted by the residency guard before an accepted result: a delayed foreign CUDA attachment was confirmed. Only the owned process was stopped; no peer was signaled. The [second interruption record](../../../../validation/gpu_exclusivity_2026-10-23/second_residency_interruption.json) preserves the exact identities. This attempt is terminal and will not be automatically reused.\n')
(new/'README.md').write_text('# Fresh continuation after a guarded residency interruption\n\nThe [previous attempt](../baseline__design_pmdi_ethyl__occ_2_3_exclusive_resume/README.md) was stopped by its residency guard before an accepted result. The peer was untouched. This fresh campaign retains the original complete evaluation 23 as an unqualified warm start, with identical chemical/numerical settings, validated controls, host allowance and final QC requirements.\n\n[The plan](plan.json) binds [the interruption proof](../../../../validation/gpu_exclusivity_2026-10-23/second_residency_interruption.json), source checkpoint and controls. The finite CPU observer remains required. Placement is restricted to an exclusively free suitable H100 with sufficient shared RAM and CPU; current device occupancy can delay an actual start. Preparation does not qualify an optimized energy, a Hessian or G.\n')
registry=base/'baseline_recoveries.json'
with registry.with_suffix('.lock').open('a+') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);x=json.loads(registry.read_text());assert x['overrides']['baseline__design_pmdi_ethyl__occ_2_3']['campaign']==str(old)
 history=base/'baseline_observer_history'/f'registry_before_second_exclusive_checkpoint_{time.time_ns()}.json';shutil.copyfile(registry,history)
 x['overrides']['baseline__design_pmdi_ethyl__occ_2_3']=dict(campaign=str(new),hosts=['lambda13'],reason='Second delayed foreign attachment detected and own worker stopped by guard. Fresh original evaluation-23 continuation, same scientific settings and host allowance; actual start only on an exclusively free suitable device.')
 x.update(previous_registry=str(history),updated=time.time());write(registry,x)
print(json.dumps(dict(status='fresh_continuation_registered',campaign=str(new),accepted_E=False,accepted_G=False,owned_predecessor_stopped=True,peer_preserved=True)))
