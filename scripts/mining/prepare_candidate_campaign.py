#!/usr/bin/env python
"""Write an unsubmitted candidate campaign with explicit preparation gates."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from redox.core.common import load_config
from redox.core.protocol import ACTIVE_SP
from redox.build.contracts import validate_group

ORDER=['exp_phthalimide_control','exp_viologen_methoxyethyl','exp_aq_1oxy','exp_ndi_hexyl_pentyloxy']


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args()
    if a.out.exists():ap.error('choose a fresh campaign directory')
    prepared=json.loads((ROOT/'results/discovery/followup_inputs_2026-10-11/batch_config.json').read_text())
    custom=load_config('custom_candidates');registered={g['id']:g for g in custom.GROUPS}
    definitions={g['id']:g for g in prepared};definitions.update(registered)
    py='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python';out=a.out.resolve()
    entries=[]
    for gid in ORDER:
        g=definitions[gid];validate_group(custom.SCAFFOLD,g)
        runs=out/'runs'/gid
        argv=[py,str(ROOT/'scripts/pipeline/run_motifs.py'),'--config',str(ROOT/'config/custom_candidates.py'),
              '--only',gid,'--out',str(runs),'--thermal-definition','native','--threads','2','--nconf','3']
        archive=[py,str(ROOT/'scripts/pipeline/stage_uma_recalculation.py'),'--only',gid,
                 '--archive-root',str(out/'archives'/gid)]
        entries.append(dict(id=gid,registered=gid in registered,states=g['states'],name=g['name'],frag=g['frag'],
            provenance=g.get('provenance'),chemistry_notes=g.get('chemistry_notes'),
            readiness='eligible_for_reviewed_UMA_recalculation' if gid in registered else 'pending_root_registration',
            prerequisites=['Use every currently free suitable GPU; preflight each node and assign disjoint CPU sets before fresh system reservations',
                           'Native thermal implementation controls pass and software versions match',
                           'Register exact reviewed input under custom_candidates.py before execution' if gid not in registered else
                           'Archive existing unverifiable legacy UMA states before recalculation; no provenance backfill'],
            archive_plan_argv=archive,archive_execute_argv=archive+['--execute'],
            bounded_uma_after_archiving=f'gpu_reserve run GPU_INDEX -- timeout --signal=TERM --kill-after=60 3h taskset -c CPU_PAIR {py} -m redox.qm.uma --only {gid}',
            pipeline_plan_argv=argv,pipeline_execute_argv=argv+['--execute','--gpu','GPU_INDEX'],
            success_conditions=['Every declared structure, charge and geometry stereochemistry passes readiness',
                'Native addressed thermal corrections and active SMD energies exist for candidate and Fc',
                'All potentials use matching candidate/Fc thermal definitions',
                'Small or disputed spin decisions receive active-protocol DFT checks',
                'New kinetic objectives remain missing/review until the physical reorganization gate passes',
                'Audited arithmetic and selected-candidate completion QC both pass'],
            inherited_parent_metrics=False))
    remaining=[definitions[gid] for gid in ORDER if gid not in registered]
    a.out.mkdir(parents=True)
    (out/'registration_proposals.json').write_text(json.dumps(remaining,indent=2)+'\n')
    plan=dict(status='prepared_not_submitted',production_order=ORDER,resource_policy='all currently free suitable GPUs within the authorized workload',
              threads_per_candidate_worker=2,active_energy_protocol=ACTIVE_SP,thermal_definition='native_solution_dft_qrrho_spin',
              quantum_jobs_submitted=0,queue=entries,
              source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'config/custom_candidates.py',ROOT/'scripts/pipeline/run_motifs.py',ROOT/'scripts/pipeline/stage_uma_recalculation.py']})
    (out/'campaign.json').write_text(json.dumps(plan,indent=2)+'\n')
    lines=['Unsubmitted campaign. GPU_INDEX and CPU_PAIR require fresh system reservations, validated GPU memory and disjoint CPU allocations. Use all currently free suitable GPUs without an arbitrary worker cap.']
    for e in entries:
        lines += ['\n'+e['id']+' — '+e['readiness'], 'Archive inspection: '+shlex.join(e['archive_plan_argv']),
                  'Pipeline plan: '+shlex.join(e['pipeline_plan_argv'])]
    (out/'command_plan.txt').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(out=str(out),status=plan['status'],ids=ORDER),indent=2))
if __name__=='__main__':main()
