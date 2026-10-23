#!/usr/bin/env python
"""Launch the prepared remaining CREST panel plus two reference-repeat searches.

Static disjoint tasks use a source-pinned worker with atomic host CPU allocation.
No GPU is used. The earlier reference pilots remain untouched; the repeats test
realized-search sensitivity for the empty backbone, which previously moved the
occupancy cycle substantially. Statistical RNG independence is not assumed.
"""
import argparse,hashlib,json,shlex,subprocess
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json

PYTHON='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python'
BASE=ROOT/'results/validation/conformer_sampling_2026-10-23'
REMOTE='''import json,os,subprocess,time
from pathlib import Path
c=json.load(__import__('sys').stdin);records=[]
for item in c['items']:
 env=dict(os.environ,PYTHONPATH=c['root']+'/src',OMP_NUM_THREADS='14',OMP_THREAD_LIMIT='14',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',XTB_THREADS='1')
 with open(item['log'],'x') as log:
  p=subprocess.Popen(item['argv'],cwd=c['root'],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 records.append(dict(item,pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],started=time.time()))
 time.sleep(.3)
print(json.dumps(dict(host=os.uname().nodename,launches=records)))
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--node',default='lambda13');args=parser.parse_args()
    preflight=json.loads((BASE/(args.node+'_preflight.json')).read_text())
    if preflight.get('ok') is not True:raise ValueError('Node preflight failed')
    dry=ROOT/'calcs/polymer_models/2026-10-23_conformer_revision/dry_lambda13_largest/status.json'
    if json.loads(dry.read_text())['status']!='dry_settings_check_complete_no_sampling':raise ValueError('Node CREST settings/dependencies not checked')
    pilot=json.loads((BASE/'pilot_start_verification.json').read_text())
    if not pilot['actual_MTD_started'] or pilot['OpenBLAS_warning_count']!=0:raise ValueError('Execution strategy has not passed startup checks')
    tasks=json.loads((BASE/'inputs/plan.json').read_text())['tasks'][1:]
    repeat=BASE/'reference_repeats';repeat.mkdir()
    for label in ('all_same','alternating'):
        old=ROOT/'results/validation/conformer_revision_2026-10-22'/label
        directory=repeat/label;directory.mkdir()
        case=json.loads((old/'case.json').read_text());case.update(mult=1,key=label+'__occ_none__parallel_repeat',
            role='Empty-backbone reference-repeat search; no statistical RNG independence claim',
            initial_geometry_source=str(old/'input.xyz'),initial_geometry_source_sha256=hashlib.sha256((old/'input.xyz').read_bytes()).hexdigest())
        atomic_write_json(directory/'case.json',case);(directory/'input.xyz').write_bytes((old/'input.xyz').read_bytes())
        tasks.append(dict(key=case['key'],input=str(directory/'input.xyz'),case=str(directory/'case.json'),
            raw_out=str(ROOT/'calcs/polymer_models/2026-10-23_conformer_revision'/case['key'])))
    if preflight['load']+14*len(tasks)>preflight['online_cpus']:raise ValueError('Node cannot cover the full requested CPU wave')
    items=[]
    for task in tasks:
        if Path(task['raw_out']).exists():raise ValueError('Preserve already started task instead of duplicate launch')
        for key in ('input','case'):task[key]=str((ROOT/task[key]).resolve())
        argv=[PYTHON,str(ROOT/'scripts/pipeline/sample_polymer_conformers.py'),'--input',task['input'],'--case',task['case'],
            '--out',task['raw_out'],'--threads','14','--max-hours','48']
        items.append(dict(key=task['key'],raw_out=task['raw_out'],argv=argv,log=str(BASE/(task['key']+'_launch.log'))))
    atomic_write_json(BASE/'remaining_wave_plan.json',dict(tasks=tasks,host=args.node,workers=len(tasks),requested_CPU_slots=14*len(tasks),
        scientific_unique_graphs_across_original_panel=12,additional_empty_reference_repeats=2,
        native_thread_caps=dict(OpenMP=14,BLAS=1),GPU_calculations=False,scope='Low-level sampling only; DFT/ensemble G pending'))
    response=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',args.node,shlex.join([PYTHON,'-c',REMOTE])],
        input=json.dumps(dict(root=str(ROOT),items=items)),text=True,capture_output=True,check=True,timeout=30)
    record=json.loads(response.stdout);atomic_write_json(BASE/('launch_'+args.node+'.json'),record)
    print(json.dumps(dict(host=args.node,launched=len(record['launches']),CPU_slots_requested=14*len(record['launches']))))


if __name__=='__main__':main()
