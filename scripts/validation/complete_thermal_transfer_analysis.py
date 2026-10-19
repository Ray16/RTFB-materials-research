#!/usr/bin/env python
"""Bounded CPU-only finalizer for preserved thermal calculations and DFT reference."""
import argparse,hashlib,json,os,socket,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT,DFT,read_result
from redox.qm.thermal import load_solution_thermal_record
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-wait-hours',type=float,default=6);a=ap.parse_args()
    p=BASE/'thermal_transfer_finalization.json'
    status=dict(status='waiting_for_models_and_reference',host=socket.gethostname(),pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().split()[21],cpu_affinity=sorted(os.sched_getaffinity(0)),
        started=time.time(),max_wait_hours=a.max_wait_hours,production_modified=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    def save():
        t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
    save()
    try:
        deadline=time.time()+a.max_wait_hours*3600
        while time.time()<deadline:
            q=BASE/'uniform_solution_thermal/status.json';h=BASE/'numerical_solution_hessian_gridresponse_verified/result.json'
            if q.exists():
                model=json.loads(q.read_text())
                if model['status']=='failed':raise RuntimeError('Uniform thermal model worker failed; preserve and inspect its exact error')
                if model['status']=='model_validation_complete' and h.exists():break
            coordinator=BASE/'numerical_solution_hessian_gridresponse_verified/coordinator_status.json'
            if coordinator.exists() and json.loads(coordinator.read_text())['status']=='failed':raise RuntimeError('Conservative DFT Hessian validation failed')
            time.sleep(15)
        else:raise RuntimeError('Thermal/reference completion deadline expired')
        reference=json.loads(h.read_text())
        if reference['status']!='numerical_reference_verified' or not reference['stationary_minimum']:
            raise RuntimeError('DFT reference is not a verified stationary minimum')
        valid=[];missing=[]
        for state in model['states']:
            raw=read_result(state['id'],state['state'],raw=True)
            rec=load_solution_thermal_record(state['selected'],DFT/state['id']/state['state']/'opt.xyz',raw['charge'],raw['mult'])
            (valid if rec is not None else missing).append(dict(id=state['id'],state=state['state'],record=state['selected']))
        d3=json.loads((BASE/'d3_correct_spin_solution_thermal_v2/status.json').read_text())
        if d3['status']!='model_validation_complete':raise RuntimeError('Current-spin d3 paired models incomplete')
        for state in d3['states']:
            raw=read_result(state['id'],state['state'],raw=True)
            if load_solution_thermal_record(state['selected'],DFT/state['id']/state['state']/'opt.xyz',raw['charge'],raw['mult']) is None:
                raise RuntimeError('Current-spin d3 thermal identity/minimum failed')
        dest=BASE/'transfer_analysis_final'
        log=BASE/'transfer_analysis_final.log'
        status.update(status='analyzing',valid_stationary_states=len(valid),invalid_model_states=missing);save()
        with log.open('w') as f:
            proc=subprocess.run([sys.executable,str(ROOT/'scripts/validation/analyze_thermal_transfer.py'),'--out',str(dest)],
                cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=300)
        if proc.returncode:raise RuntimeError('Transfer analysis failed; inspect exact log')
        transfer=json.loads((dest/'transfer_reference.json').read_text())
        if not transfer['reference_available']:raise RuntimeError('Transfer analysis omitted verified reference')
        # Report calibrated measurements; never invent general validation from one state.
        status.update(status='transfer_analysis_completed' if not missing else 'transfer_analysis_completed_with_model_failures',
            output=str(dest),reference=str(h),valid_model_states=valid,model_failures=missing,
            single_state_transfer_calibration=transfer['model_error_at_aq_red2_eV'],
            general_transfer_accuracy_validated=False,finished=time.time())
        save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
