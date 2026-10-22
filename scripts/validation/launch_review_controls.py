#!/usr/bin/env python
"""Launch the finite, preflighted review controls through the shared GPU gate."""
import json,shlex,subprocess
from pathlib import Path
from redox.core.common import ROOT

PYTHON='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
BASIS=ROOT/'results/validation/neutral_basis_2026-10-22'
FINAL=ROOT/'results/validation/saved_trimer_final_2026-10-22'
REMOTE='''import json,os,subprocess,sys,time
from pathlib import Path
cfg=json.load(sys.stdin);root=Path(cfg['root']);result=[]
for item in cfg['items']:
 log=Path(item['log']);env=dict(os.environ,PYTHONPATH=str(root/'src'))
 for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
 with log.open('x') as stream:
  child=subprocess.Popen(item['argv'],cwd=root,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
 result.append(dict(item,pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],started=time.time()))
print(json.dumps(dict(host=os.uname().nodename,launches=result)))
'''


def main():
    plan=json.loads((BASIS/'plan.json').read_text())
    # All five allowed controls, across all six suitable slots observed for this
    # six-task wave. GPU reservations recheck occupancy atomically at start.
    assignments=[('lambda13',0,'saved_trimer'),('lambda13',1,plan['tasks'][0]['key']),
        ('lambda13',2,plan['tasks'][1]['key']),('lambda13',3,plan['tasks'][2]['key']),
        ('lambda11',0,plan['tasks'][3]['key']),('lambda11',1,plan['tasks'][4]['key'])]
    for host in ('lambda13','lambda11'):
        preflight=json.loads((BASIS/(host+'_preflight.json')).read_text())
        if preflight.get('ok') is not True:raise ValueError('Node preflight failed')
        items=[]
        for node,gpu,key in assignments:
            if node!=host:continue
            folder=FINAL if key=='saved_trimer' else BASIS
            destination=folder/'execution'/key
            if destination.exists():raise ValueError('Preserve previous execution rather than duplicate launch')
            command=([PYTHON,str(ROOT/'scripts/validation/finish_saved_trimer.py'),'--plan',str(FINAL/'plan.json')]
                     if key=='saved_trimer' else [PYTHON,str(ROOT/'scripts/validation/neutral_basis_sensitivity.py'),'--campaign',str(BASIS),'--task',key])
            argv=[PYTHON,'-m','redox.workflows.gpu','--gpu',str(gpu),'--out',str(destination),'--',*command]
            items.append(dict(key=key,gpu=gpu,argv=argv,log=str(folder/(key+'_launch.log'))))
        argv=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join([PYTHON,'-c',REMOTE])]
        response=subprocess.run(argv,input=json.dumps(dict(root=str(ROOT),items=items)),text=True,capture_output=True,timeout=30,check=True)
        result=json.loads(response.stdout);(BASIS/('launch_'+host+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(dict(host=host,launched_wrappers=len(result['launches']))))


if __name__=='__main__':main()
