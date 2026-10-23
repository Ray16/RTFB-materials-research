#!/usr/bin/env python
"""Publish source-current candidate evidence and feasibility-first decision views.

Read-only scientific reconstruction and figure generation only: no QM submission,
new candidate registration, parent-property inheritance or missing-value fill.
Stages use separate interpreters so CPU affinity and native libraries do not leak.
"""
import argparse,hashlib,importlib.util,json,os,shutil,subprocess,sys,time
from datetime import datetime,timezone
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(argv,*,allocate=False,threads=1):
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]=str(threads)
    if allocate:
        # The stage itself does not allocate resources. Bind before importing it;
        # do not use this wrapper for stages that own their own allocation.
        code='import os,sys; from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity; '\
             'guard=allocation_guard(); guard.__enter__(); host_capacity(1); os.sched_setaffinity(0,cpu_allocation(1)); guard.__exit__(None,None,None); os.execv(sys.executable,[sys.executable,*sys.argv[1:]])'
        argv=[argv[0],'-c',code,*argv[1:]]
    subprocess.run(argv,cwd=ROOT,env=env,check=True)


def current_benchmark():
    from redox.validation.benchmark_panel import review_report
    panel_path=ROOT/'results/validation/publication_benchmark/panel_v1.json';panel=json.loads(panel_path.read_text())
    current=ROOT/'results/validation/publication_benchmark/current.json';previous=json.loads(current.read_text())
    try:
        old=next(Path(p).parent for p in previous['review_snapshot_sha256'] if Path(p).name=='summary.json')
        report=review_report(panel,old)
    except (OSError,ValueError,KeyError,StopIteration):
        stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
        review=ROOT/'results/reviews'/datetime.now(timezone.utc).strftime('%Y-%m-%d')/('publication_pipeline_'+stamp)
        ids=','.join(sorted({r['id'] for r in panel['targets']}))
        run([sys.executable,str(ROOT/'scripts/validation/benchmark/review_matched_parents.py'),'--only',ids,'--out',str(review)],threads=2)
        report=review_report(panel,review)
    if report['admitted_events']!=report['target_events']:raise ValueError('Fixed publication benchmark is incomplete')
    atomic_write_json(current,report)
    old_values={(r['id'],r['event'],r['doi']):r['calc_V'] for r in previous['outliers']}
    new_values={(r['id'],r['event'],r['doi']):r['calc_V'] for r in report['outliers']}
    changed=old_values.keys()!=new_values.keys() or any(abs(v-new_values[k])>1e-9 for k,v in old_values.items())
    return dict(events=report['admitted_events'],numerical_benchmark_changed=changed,manuscript_refresh_required=changed)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--feed',type=Path,default=ROOT/'results/candidates/current.json')
    parser.add_argument('--out',type=Path);parser.add_argument('--plot-python',type=Path,default=ROOT/'.venv-figures/bin/python')
    parser.add_argument('--comparison-already-updated',action='store_true');args=parser.parse_args()
    if not args.plot_python.is_file():raise ValueError('A plotting interpreter must be supplied')
    benchmark=current_benchmark()
    helper=ROOT/'scripts/plotting/candidates/verified_overlay.py'
    spec=importlib.util.spec_from_file_location('decision_feed_freshness',helper);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    feed=json.loads(args.feed.read_text())
    if not module.snapshot_inputs_current(feed['summary']):
        canonical=ROOT/'results/candidates/current.json'
        if args.feed.resolve()!=canonical.resolve():raise ValueError('Explicit feed is stale; preserve it and collect a new snapshot')
        run([sys.executable,str(ROOT/'scripts/pipeline/collect_verified_candidates.py'),'--out',str(args.feed.resolve().parent),'--max-hours','1'],allocate=True)
        feed=json.loads(args.feed.read_text())
    if not module.snapshot_inputs_current(feed['summary']):raise ValueError('Candidate evidence is not source-current')
    before=sha(args.feed)
    out=(args.out or ROOT/'results/selection/campaigns'/datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')).resolve()
    if not out.is_relative_to(ROOT):raise ValueError('Publish result metadata inside the project; portable exports are separate')
    run([sys.executable,str(ROOT/'scripts/pipeline/collect_candidate_decisions.py'),'--feed',str(args.feed),'--out',str(out)])
    if not args.comparison_already_updated:
        run([str(args.plot_python),str(ROOT/'scripts/plotting/candidates/update_comparison_collection.py'),'--feed',str(args.feed)],allocate=True)
    else:
        manifest=json.loads((ROOT/'results/candidates/figures/metadata/manifest.json').read_text())
        if manifest['snapshot']!=feed['snapshot_path']:raise ValueError('Comparison figures do not match the declared current feed')
    for script in ('plot_descriptor_pareto.py','plot_candidate_axis_matrix.py'):
        run([str(args.plot_python),str(ROOT/'scripts/plotting/candidates'/script),'--decisions',str(out)])
    if sha(args.feed)!=before or not module.snapshot_inputs_current(feed['summary']):raise ValueError('Candidate evidence changed during publication')
    summary=json.loads((out/'summary.json').read_text())
    source=summary['source_sha256']
    if any(sha(p)!=h for p,h in source.items()):raise ValueError('Decision inputs changed during publication')
    figures=out/'figure_snapshot';figures.mkdir()
    for stem in ('descriptor_pareto','candidate_axis_matrix'):
        for extension in ('pdf','png','svg','json'):
            path=ROOT/'results/figures/selection'/(stem+'.'+extension)
            shutil.copyfile(path,figures/path.name)
    for name in ('plot_descriptor_pareto.py','plot_candidate_axis_matrix.py','candidate_highlights.py','verified_overlay.py'):
        path=ROOT/'scripts/plotting/candidates'/name;target=out/'source_snapshot/producers'/path.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(path.read_bytes())
    target=out/'source_snapshot/producers/scripts/pipeline/publish_candidate_decisions.py'
    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(Path(__file__).read_bytes())
    outputs={str(p.relative_to(ROOT)):sha(p) for p in out.rglob('*') if p.is_file()}
    atomic_write_json(out/'publication_manifest.json',dict(kind='source_current_all_axes_publication',created=time.time(),
        source_sha256=source,output_sha256=outputs,benchmark=benchmark,quantum_jobs_submitted=0,
        publisher_sha256=sha(__file__),scope='Decision criteria/evidence may remain incomplete; eligible means conditional on the declared profile.'))
    pointer=dict(guide='results/selection/README.md',report=str(out.relative_to(ROOT)),summary=str((out/'summary.json').relative_to(ROOT)),
        decisions=str((out/'decisions.json').relative_to(ROOT)),fronts=str((out/'decision_fronts.json').relative_to(ROOT)),
        manifest=str((out/'publication_manifest.json').relative_to(ROOT)),status='criteria_and_evidence_incomplete' if not summary['eligible_ids'] else 'conditional_shortlist_published',
        eligible=len(summary['eligible_ids']),input_candidate_as_of=feed['as_of_Chicago'])
    atomic_write_json(ROOT/'results/selection/current.json',pointer)
    print(json.dumps(pointer))


if __name__=='__main__':main()
