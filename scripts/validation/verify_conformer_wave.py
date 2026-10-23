#!/usr/bin/env python
"""Read-only verification of actual CPU CREST worker/child ownership and progress."""
import argparse,hashlib,json,shlex,subprocess,time
from pathlib import Path
from redox.core.common import ROOT

REMOTE='''import json,os
from pathlib import Path
items=json.load(__import__('sys').stdin);rows=[];used=set()
for item in items:
 p=Path(item['raw_out']);state=json.loads((p/'status.json').read_text());processes=[]
 for key in ('pid','predictor_pid'):
  pid=state[key];proc=Path(f'/proc/{pid}')
  if not proc.exists():continue
  start=(proc/'stat').read_text().rsplit(')',1)[1].split()[19]
  expected=state['process_start'] if key=='pid' else state['predictor_start']
  if start!=expected:continue  # A reused PID is not the old worker/predictor.
  caps={}
  for field in (proc/'environ').read_bytes().split(bytes([0])):
   if b'=' not in field:continue
   k,v=field.split(b'=',1)
   if k.decode() in state['thread_caps']:caps[k.decode()]=v.decode()
  affinity=sorted(os.sched_getaffinity(pid))
  if affinity!=state['CPU_affinity'] or caps!=state['thread_caps']:raise ValueError('Affinity/thread caps differ')
  processes.append(dict(role=key,pid=pid,process_start=start,CPU_affinity=affinity,thread_caps=caps))
 if processes and used.intersection(state['CPU_affinity']):raise ValueError('Workers overlap CPU sets')
 if processes:used.update(state['CPU_affinity'])
 log=(p/'crest.log').read_text()
 roles={r['role'] for r in processes}
 lifecycle=('owner_and_predictor_live' if roles=={'pid','predictor_pid'} else 'owner_live_predictor_terminal_qualification_pending' if 'pid' in roles else 'owner_absent_predictor_live_preserve_work' if 'predictor_pid' in roles else 'no_live_handle_inspect_terminal_evidence')
 rows.append(dict(key=item['key'],status=state['status'],lifecycle=lifecycle,processes=processes,actual_MTD_started='starting MTD' in log,
  trajectory_files=len(list((p/'MDFILES').glob('*.trj'))),log_bytes=(p/'crest.log').stat().st_size,
  OpenBLAS_warning_count=log.count('OpenBLAS Warning'),DFT_reranking_complete=state['DFT_reranking_complete'],ensemble_G_computed=state['ensemble_G_computed']))
print(json.dumps(dict(host=os.uname().nodename,workers=rows,disjoint_CPU_count=len(used))))
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--launch',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    launch=json.loads(args.launch.read_text());host=launch['host'];items=[dict(key=r['key'],raw_out=r['raw_out']) for r in launch['launches']]
    python='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
    result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join([python,'-c',REMOTE])],
        input=json.dumps(items),capture_output=True,text=True,check=True,timeout=30)
    proof=json.loads(result.stdout);proof.update(checked=time.time(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Actual sampler starts/resource inheritance; not completed sampling, DFT or ensemble free energies')
    with args.out.open('x') as stream:json.dump(proof,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(host=host,workers=len(proof['workers']),MTD_started=sum(r['actual_MTD_started'] for r in proof['workers']),disjoint_CPUs=proof['disjoint_CPU_count'])))


if __name__=='__main__':main()
