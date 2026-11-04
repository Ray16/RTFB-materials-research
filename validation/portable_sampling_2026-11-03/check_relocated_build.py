"""Run the full manuscript validator in a relocated source tree with data paths denied."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,required=True)
parser.add_argument('--prefix',type=Path,required=True)
args=parser.parse_args();source=args.source.resolve()
outputs={suffix:Path(str(args.prefix)+suffix) for suffix in ('.stdout.log','.stderr.log','.json')}
if any(p.exists() for p in outputs.values()):raise FileExistsError('Preserve earlier build attempts')
blocked=[ROOT,ROOT.parent/'archive',*[ROOT.parent/name for name in
    ('redox_candidate_continuations','redox_kernel_controls','redox_polymer_continuations',
     'redox_counterion_controls','redox_motif_interactions')]]
code='''import pathlib,runpy,sys
root=pathlib.Path(sys.argv[1]);blocked=[pathlib.Path(x).resolve() for x in sys.argv[2:]]
sys.path[:0]=[str(root/'manuscript'),str(root/'src')]
def deny(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes)):
  path=pathlib.Path(args[0]).resolve()
  if any(path==b or path.is_relative_to(b) for b in blocked):raise PermissionError('Original scientific path blocked: '+str(path))
sys.addaudithook(deny)
runpy.run_path(str(root/'manuscript/build.py'),run_name='__main__')
'''
# The builder allocates its own CPU; the waiting parent must not reserve that CPU.
run=subprocess.run([sys.executable,'-I','-c',code,str(source),*map(str,blocked)],
                   cwd=source,capture_output=True,text=True)
outputs['.stdout.log'].write_text(run.stdout);outputs['.stderr.log'].write_text(run.stderr)
report=dict(status='relocated_full_manuscript_build_passed' if run.returncode==0 else 'relocated_build_failed_preserved',
    returncode=run.returncode,source=str(source),blocked_roots=list(map(str,blocked)),
    original_scientific_paths_blocked=True,new_quantum_validation=False)
if run.returncode==0:
    built=json.loads((source/'manuscript/build_check.json').read_text())
    report['documents']={name:row['pages'] for name,row in built['documents'].items()}
    report['occupied_sampling_check']=built['occupied_sampling_check']
    report['ionic_preparation_check']=built['ionic_preparation_check']
outputs['.json'].write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],returncode=run.returncode,error_tail=run.stderr[-1600:])))
raise SystemExit(run.returncode)
