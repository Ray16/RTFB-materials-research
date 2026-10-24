"""Stop only this exact redox worker after verifying a shared-device conflict."""
from pathlib import Path
import json,os,signal,socket,subprocess,time
import psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');C=R/'results/polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume'
row=json.loads((C/'energy/workers/lambda13_gpu3_queue.json').read_text());p=psutil.Process(row['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19];assert start==row['process_start'];assert str(C) in ' '.join(p.cmdline());assert p.pid==1925616
raw=R/'calcs/polymer_models/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume/energy/baseline__design_pmdi_ethyl__occ_2_3__pose0/rest/geometry_checkpoints'
latest=json.loads((raw/'latest.json').read_text());assert latest['SCF_converged'];assert (raw/'evaluated_trajectory.xyz').stat().st_size>0
lines=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_memory','--format=csv,noheader,nounits'],text=True).splitlines();rows=[r.strip().split(', ') for r in lines if r.strip()];ours=[r for r in rows if int(r[0])==p.pid];assert len(ours)==1;others=[r for r in rows if r[1]==ours[0][1] and int(r[0])!=p.pid];assert others
proof=dict(status='stopping_owned_worker_for_device_exclusivity',checked=time.time(),host=socket.gethostname(),pid=p.pid,process_start=start,GPU=row['gpu'],resident_processes_on_device=[r for r in rows if r[1]==ours[0][1]],latest_saved_unqualified_evaluation=latest,unrelated_processes_signaled=False)
out=R/'results/validation/pcm_first_response_workspace_v2_2026-10-23/device_exclusivity_recovery.json';out.write_text(json.dumps(proof,indent=2)+'\n')
p.send_signal(signal.SIGTERM)
for _ in range(200):
 if not p.is_running() or p.status()==psutil.STATUS_ZOMBIE:break
 time.sleep(.1)
else:raise RuntimeError('Owned worker did not stop; preserve claims and inspect before another signal')
proof.update(status='owned_worker_stopped_checkpoints_preserved',stopped=time.time(),worker_alive=False)
for r in others:
 assert psutil.pid_exists(int(r[0]))
proof['unrelated_processes_still_present']=True;out.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
