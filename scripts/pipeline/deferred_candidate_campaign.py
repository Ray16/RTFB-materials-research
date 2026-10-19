#!/usr/bin/env python
"""Prepare or run one deferred candidate worker after GPU1's physical slot releases.

Default only writes a concrete queued manifest. --execute waits at most12h, performs a
fresh Lambda preflight, and executes each candidate sequentially with2CPU threads and
atomic reservation through run_motifs.py. No legacy UMA provenance is invented.
"""
import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import shutil
import socket
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from redox.core.common import load_config,read_manifest,read_result,free_energy,LIBRARY,UMA
from redox.build.contracts import validate_group
GATE=Path('/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve')
POLICY=Path('/nfs/lambda_stor_01/homes/rzhu/.config/gpu_reserve/config.env')
PREFLIGHT=Path('/nfs/lambda_stor_01/homes/rzhu/.codex-personal/skills/lambda-fleet/scripts/preflight_node.py')
DEFAULT_IDS=['exp_phthalimide_control','exp_viologen_methoxyethyl','exp_aq_1oxy','exp_ndi_hexyl_pentyloxy']
PY=Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python')


def start_identity(pid):
    return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]


def save(path,value):
    staged=path.with_name(path.name+f'.{os.getpid()}.tmp');staged.write_text(json.dumps(value,indent=2)+'\n');staged.replace(path)


def run_child(argv,log,env,seconds):
    with log.open('w') as stream:
        child=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            return child.wait(timeout=seconds)
        except BaseException:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
            raise


def cpu_slots(n):
    # A controller pinned to one CPU cannot allocate a two-CPU child from its own
    # narrowed affinity. This brief single-thread probe expands ONLY its own
    # affinity to the host/cgroup-allowed set and invokes the shared allocator.
    path=ROOT/'scripts/pipeline/run_motifs.py'
    code=("import os,json,importlib.util; os.sched_setaffinity(0,range(os.cpu_count())); "
          f"s=importlib.util.spec_from_file_location('allocation',{str(path)!r}); "
          "m=importlib.util.module_from_spec(s);s.loader.exec_module(m); "
          f"print(json.dumps(m.cpu_allocation({n})))")
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    try:
        return json.loads(subprocess.check_output([str(PY),'-c',code],env=env,text=True,timeout=30,stderr=subprocess.PIPE))
    except (subprocess.CalledProcessError,subprocess.TimeoutExpired,ValueError) as exc:
        raise RuntimeError('Fresh host CPU allocation failed') from exc


def native_fc_pair_qc():
    reference=[]
    for st in ('neu','ox'):
        r=read_result('ferrocene',st) or {};G,qc=free_energy(r)
        reference.append(dict(state=st,thermal_definition=r.get('thermal_definition'),thermal_qc=qc,
            thermal_protocol_hash=r.get('thermal_protocol_hash'),active_sp=r.get('sp_status'),free_energy_present=G is not None))
    ready=all(r['thermal_definition']=='native_solution_dft_qrrho_spin' and r['thermal_qc']=='ok'
              and r['thermal_protocol_hash'] and r['active_sp']=='ok' and r['free_energy_present'] for r in reference)
    return bool(ready),reference


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--ids',default=','.join(DEFAULT_IDS));ap.add_argument('--execute',action='store_true')
    ap.add_argument('--wait-seconds',type=int,default=43200);ap.add_argument('--candidate-timeout-seconds',type=int,default=43200)
    a=ap.parse_args();ids=a.ids.split(',')
    if not 0<a.wait_seconds<=43200 or not 0<a.candidate_timeout_seconds<=43200:ap.error('deadlines must be positive and at most12h')
    if len(ids)!=len(set(ids)):ap.error('duplicate candidate IDs')
    cfg=load_config('custom_candidates');byid={g['id']:g for g in cfg.GROUPS}
    if set(ids)-set(byid):ap.error('every selected input must already be registered')
    for gid in ids:validate_group(cfg.SCAFFOLD,byid[gid])
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    status_path=out/'status.json'
    planned=dict(status='queued_not_running',ids=ids,expected_states=sum(len(byid[g]['states']) for g in ids),
       host='lambda1',gpu=1,max_total_gpu_workers=2,max_candidate_workers=1,threads=2,
       controller_cpu=6,controller_threads=1,
       wait_seconds=a.wait_seconds,candidate_timeout_seconds=a.candidate_timeout_seconds,
       registration=[byid[g] for g in ids],queue=[],controller_script=str(Path(__file__).resolve()),
       interpreter=str(PY),policy=str(POLICY),created=time.time(),quantum_jobs_submitted=0,
       quantum_execution_status='not_started',
       source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
          [Path(__file__),ROOT/'config/custom_candidates.py',ROOT/'scripts/pipeline/run_motifs.py',ROOT/'scripts/pipeline/stage_uma_recalculation.py']})
    planned['quantum_jobs_submitted_at_preparation']=planned.pop('quantum_jobs_submitted')
    if not a.execute:
        if status_path.exists():ap.error('queue manifest already exists; use --execute or choose a new output')
        save(status_path,planned);print(status_path);return
    if socket.gethostname().split('.')[0]!='lambda1':ap.error('this reviewed controller is assigned to lambda1')
    if os.sched_getaffinity(0)!={6}:ap.error('launch controller with taskset -c6; children receive their own fresh allocation')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
    policy=POLICY.read_text()
    import shlex
    parsed={k:shlex.split(v)[0] for line in policy.splitlines() if line.strip() and not line.lstrip().startswith('#') and '=' in line for k,v in [line.split('=',1)]}
    if 'lambda1' not in parsed.get('ALLOWED_HOSTS','').split() or 'lambda1' in parsed.get('RESERVED_HOSTS','').split():ap.error('current reservation policy excludes lambda1')
    lock=(out/'.controller.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    status=json.loads(status_path.read_text()) if status_path.exists() else planned
    if status.get('ids')!=ids:ap.error('queue identities changed; choose a new run output')
    for path,digest in status['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:ap.error('reviewed source changed; regenerate queue manifest before execute')
    status.update(status='waiting_for_physical_slot',pid=os.getpid(),process_start=start_identity(os.getpid()),started=time.time());save(status_path,status)
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'));env['PATH']=str(PY.parent)+os.pathsep+env.get('PATH','')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='2'
    import site
    paths=[str(p) for base in site.getsitepackages() for p in Path(base).glob('nvidia/*/lib')]
    env['LD_LIBRARY_PATH']=':'.join(paths+[env.get('LD_LIBRARY_PATH','')])
    child=None;claim=None
    def stop(signum,frame):raise KeyboardInterrupt(f'signal{signum}')
    signal.signal(signal.SIGTERM,stop)
    try:
        deadline=time.time()+a.wait_seconds
        while time.time()<deadline:
            free=subprocess.check_output([str(GATE),'list','--hosts','lambda1'],text=True)
            if 'lambda1 1' in free.splitlines():
                try:cpus=cpu_slots(2)
                except RuntimeError as exc:
                    status.update(status='waiting_for_cpu_allocation',resource_issue=str(exc));save(status_path,status)
                    time.sleep(15);continue
                break
            time.sleep(15)
        else:raise RuntimeError('physical slot/resource deadline expired')
        preflight=[str(PY),str(PREFLIGHT),'lambda1','--python',str(PY)]
        for module in ('rdkit','ase','pyscf','gpu4pyscf','fairchem.core','redox.qm.native_thermal'):preflight+=['--module',module]
        preflight+=['--path',str(ROOT),'--bin',str(GATE)]
        xtb=env.get('XTB_BIN') or shutil.which('xtb',path=env['PATH'])
        if not xtb or not Path(xtb).is_file():raise RuntimeError('xTB executable is unavailable')
        env['XTB_BIN']=xtb;preflight+=['--bin',xtb]
        preflight+=['--library-path',env['LD_LIBRARY_PATH']]
        if run_child(preflight,out/'preflight.log',env,120):raise RuntimeError('Lambda dependency preflight failed; no candidate claimed')
        from redox.qm.uma import model_identity,DEFAULT_MODEL,seed_identity,cache_errors
        status['uma_model_identity']=model_identity(DEFAULT_MODEL)
        status['seed_identity']={gid:seed_identity(sorted((LIBRARY/gid/'conformers').glob('conf_*.xyz'))) for gid in ids}
        if any(not seeds for seeds in status['seed_identity'].values()):raise RuntimeError('Registered motif lacks shared conformer seeds')
        status['xtb_executable']=xtb
        from redox.qm.native_thermal import implementation_controls
        while True:
            try:
                status['native_implementation_proof']=implementation_controls()
                break
            except (FileNotFoundError, ValueError, RuntimeError) as exc:
                if 'version changed' in str(exc):raise
                status.update(status='waiting_for_native_controls',control_issue=str(exc));save(status_path,status)
                if time.time()>deadline:raise RuntimeError('Native implementation controls did not pass within queue deadline')
                time.sleep(15)
        while True:
            ready,reference=native_fc_pair_qc()
            status['fc_pair_qc']=reference
            if ready:break
            status['status']='waiting_for_native_Fc_pair';save(status_path,status)
            if time.time()>deadline:raise RuntimeError('Both Fc native reference states did not pass within queue deadline')
            time.sleep(15)
        save(status_path,status)
        for gid in ids:
            completed=next((e for e in reversed(status['queue']) if e['id']==gid and e.get('status')=='selected_outputs_verified'),None)
            if completed:
                from redox.validation.readiness import inspect_states
                current=inspect_states([r for r in read_manifest() if r['id']==gid])
                if current and all(r['status']=='ready' for r in current):continue
                raise RuntimeError('Previously completed state changed; inspect before a new attempt')
            cpus=cpu_slots(2)  # allocation snapshot is refreshed before each item
            if os.getloadavg()[0]+2>(os.cpu_count() or 1):raise RuntimeError('insufficient current CPU capacity')
            claim=out/f'{gid}.claim';claim.mkdir()
            save(claim/'owner.json',dict(host='lambda1',pid=os.getpid(),process_start=start_identity(os.getpid())))
            attempt=len([e for e in status['queue'] if e['id']==gid])+1
            entry=dict(id=gid,attempt=attempt,status='preparing',cpus=cpus,gpu=1,started=time.time());status['queue'].append(entry);save(status_path,status)
            archive=out/'archives'/f'{gid}_attempt{attempt}'
            argv=[str(PY),str(ROOT/'scripts/pipeline/stage_uma_recalculation.py'),'--only',gid,'--archive-root',str(archive),'--execute']
            if run_child(argv,out/f'{gid}_attempt{attempt}_archive.log',env,120):raise RuntimeError('UMA archive failed; old outputs preserved')
            cpus=cpu_slots(2);entry['cpus']=cpus
            uma_argv=[str(GATE),'run','1','--label','redox-candidate-uma-'+gid,'--','taskset','-c',','.join(map(str,cpus)),
                      'timeout','--signal=TERM','--kill-after=60','3h',str(PY),'-u','-m','redox.qm.uma','--only',gid]
            uma_log=out/f'{gid}_attempt{attempt}_uma.log'
            entry.update(status='UMA_starting',uma_argv=uma_argv)
            status.update(status='candidate_UMA_starting',quantum_execution_status='UMA_stage_requested')
            save(status_path,status)
            while True:
                rc=run_child(uma_argv,uma_log,env,10900)
                if rc==0:break
                if 'BUSY/reserved -> not running' not in uma_log.read_text():
                    raise RuntimeError('Bounded UMA recalculation failed; archive and partial outputs preserved')
                entry['status']='waiting_for_reservation';status['status']='waiting_for_reservation';save(status_path,status)
                if time.time()>deadline:raise RuntimeError('reservation contention deadline expired')
                time.sleep(15)
            verified=[]
            for row in read_manifest():
                if row['id']!=gid:continue
                sd=UMA/gid/row['state'];rec=json.loads((sd/'result.json').read_text())
                errors=cache_errors(rec,row,sd/'relaxed.xyz',sorted((LIBRARY/gid/'conformers').glob('conf_*.xyz')),
                                    status['uma_model_identity'],.02)
                if errors:raise RuntimeError('UMA output verification failed: '+','.join(errors))
                verified.append(row['state'])
            if len(verified)!=len(byid[gid]['states']):raise RuntimeError('UMA declared state coverage is incomplete')
            entry['verified_UMA_states']=verified
            status['quantum_execution_status']='UMA_state_outputs_verified';save(status_path,status)
            runout=out/'runs'/f'{gid}_attempt{attempt}'
            argv=['taskset','-c',','.join(map(str,cpus)),str(PY),str(ROOT/'scripts/pipeline/run_motifs.py'),'--config',str(ROOT/'config/custom_candidates.py'),
                  '--only',gid,'--out',str(runout),'--execute','--gpu','1','--threads','2','--cpus',','.join(map(str,cpus)),
                  '--nconf','3','--thermal-definition','native']
            entry.update(status='running',argv=argv,run_output=str(runout))
            status.update(status='candidate_running',quantum_execution_status='native_property_pipeline_requested');save(status_path,status)
            pipeline_log=out/f'{gid}_attempt{attempt}_pipeline.log'
            while True:
                rc=run_child(argv,pipeline_log,env,a.candidate_timeout_seconds)
                if rc==0:break
                mp=runout/'run_manifest.json'
                interrupted=json.loads(mp.read_text()) if mp.exists() else {}
                steps=interrupted.get('steps',[])
                gated_log=Path(steps[-1]['log']) if steps else None
                if gated_log is None or not gated_log.exists() or 'BUSY/reserved -> not running' not in gated_log.read_text():break
                if time.time()>deadline:raise RuntimeError('reservation contention deadline expired')
                # Resource contention does not consume a scientific attempt; runner
                # resumes its durable per-stage outputs into separate attempt logs.
                if '--resume' not in argv:argv.append('--resume')
                time.sleep(15)
            entry['returncode']=rc
            if rc:raise RuntimeError('Candidate pipeline failed; partial outputs preserved, inspect before resume')
            manifest=json.loads((runout/'run_manifest.json').read_text())
            if manifest.get('status') not in ('selected_outputs_verified','selected_properties_audited','tables_audited'):
                raise RuntimeError('Pipeline completion has no declared verification stage')
            from redox.validation.readiness import inspect_states
            readiness=inspect_states([r for r in read_manifest() if r['id']==gid])
            if not readiness or any(r['status']!='ready' for r in readiness):raise RuntimeError('Declared candidate state QC is incomplete after pipeline exit')
            import csv
            table_root=Path(manifest.get('attempt_output',runout))/'tables'
            score_rows=[r for r in csv.DictReader(open(table_root/'scorecard.csv')) if r['id']==gid]
            ranking_ready=bool(score_rows) and all(r.get('lambda_qc')=='ok' and r.get('status')=='candidate'
                                                  and r.get('spin_confidence') not in ('wrong_ground_state','DFT_check_required') for r in score_rows)
            entry.update(status='selected_outputs_verified',finished=time.time(),scientific_accuracy_validated=False,
                         property_ranking_ready=ranking_ready,readiness=readiness)
            save(status_path,status);shutil_claim=claim;claim=None
            shutil.rmtree(shutil_claim)
        status.update(status='candidate_sequence_verified',finished=time.time(),
                      quantum_execution_status='declared_state_outputs_verified');save(status_path,status)
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='failed_preserved_partial_outputs',error=repr(exc),finished=time.time());save(status_path,status);raise
    finally:
        if claim is not None:
            shutil.rmtree(claim)
        lock.close()
if __name__=='__main__':main()
