#!/usr/bin/env python3
"""Capture actual ionic optimizer ownership, residency and progress at one time."""
import argparse
import csv
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def checks(row,process,residents):
    owner=row['stage']['owner'];launch=row['launch'];pid=owner['pid'];caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    gpu=[p for p in residents if p['uuid']==launch['device']['uuid']]
    return dict(actual_process_live=process.get('process_start')==owner['process_start'],
        worker_command_bound=row['worker'] in process.get('argv',[]) and row['plan'] in process.get('argv',[]),
        CPU_affinity_bound=process.get('cpus')==launch['CPU_affinity'] and len(launch['CPU_affinity'])==2,
        thread_pools_capped=all(process.get('env',{}).get(k)=='2' for k in caps),
        NumPy_import_profile=process.get('env',{}).get('NUMPY_MADVISE_HUGEPAGE')=='0',
        exclusive_actual_GPU_residency=len(gpu)==1 and gpu[0]['pid']==pid,
        actual_SCF_started=row.get('SCF_progress',{}).get('cycle',0)>0,
        scientific_worker_running=row['stage']['status']=='charged_polymer_solution_optimization_running')


def observe(manifest):
    data=json.loads(manifest.read_text());rows=[];by_host={}
    for entry in data['states']:
        namespace=Path(entry['namespace'])
        if sha256_file(namespace)!=entry['namespace_sha256']:raise ValueError('Changed roster namespace')
        spec=verify(namespace);case=Path(spec['case']);stage=json.loads((case/'status.json').read_text())
        placements=list(case.parent.glob('placement*/launch.json'))
        if len(placements)!=1:raise ValueError('One actual placement required')
        launch=json.loads(placements[0].read_text());progress=case/'SCF_progress.json'
        row=dict(key=entry['key'],namespace=str(namespace),namespace_sha256=entry['namespace_sha256'],
            stage=stage,launch=launch,worker=str(Path(spec['sources'])/'worker.py'),plan=str(case/'plan.json'),
            SCF_progress=json.loads(progress.read_text()) if progress.exists() else {},
            saved_geometry_evaluations=len(list((case/'geometry_evaluations').glob('*.json'))))
        rows.append(row);by_host.setdefault(stage['owner']['host'],[]).append(stage['owner']['pid'])
    observations={}
    code='''import json,os,pathlib,subprocess,sys
result={}
for pid in json.loads(sys.argv[1]):
 p=pathlib.Path('/proc')/str(pid)
 try:
  env=dict(x.split('=',1) for x in (p/'environ').read_text().split('\\0') if '=' in x)
  result[str(pid)]=dict(process_start=(p/'stat').read_text().rsplit(')',1)[1].split()[19],cpus=sorted(os.sched_getaffinity(pid)),argv=(p/'cmdline').read_text().split('\\0'),env={k:env.get(k) for k in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS','NUMPY_MADVISE_HUGEPAGE','CUDA_VISIBLE_DEVICES']})
 except FileNotFoundError:result[str(pid)]={}
resident=subprocess.run(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_memory','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=20,check=True)
print(json.dumps(dict(processes=result,resident_csv=resident.stdout,load_1m=os.getloadavg()[0],online_cpus=os.cpu_count())))
'''
    for host,pids in by_host.items():
        argv=['python3','-c',code,json.dumps(pids)]
        run=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(argv)],capture_output=True,text=True,timeout=35,check=True)
        record=json.loads(run.stdout);record['residents']=[dict(pid=int(p),uuid=u.strip(),used_memory_MiB=int(m)) for p,u,m in csv.reader(record.pop('resident_csv').splitlines(),skipinitialspace=True)]
        observations[host]=record
    occupied={}
    for row in rows:
        host=row['stage']['owner']['host'];actual=observations[host];pid=str(row['stage']['owner']['pid'])
        row['checks']=checks(row,actual['processes'][pid],actual['residents'])
        for cpu in row['launch']['CPU_affinity']:
            if (host,cpu) in occupied:raise ValueError('Overlapping CPU allocations in the geometry panel')
            occupied[host,cpu]=row['key']
    count=sum(all(row['checks'].values()) for row in rows)
    return dict(status='ten_ionic_optimizers_actual_starts_verified' if count==10 else 'ionic_optimizer_start_observation_incomplete',
        checked=time.time(),manifest=str(manifest),manifest_sha256=sha256_file(manifest),target_states=10,verified_running_starts=count,
        states=rows,node_observations=observations,completed_geometries=0,native_G_qualified=False,material_admission=False,
        producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve previous actual-process observations')
    result=observe(a.manifest.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','verified_running_starts','target_states')}))
