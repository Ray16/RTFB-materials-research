#!/usr/bin/env python
"""Plan or explicitly archive unverifiable UMA caches before a bounded recalculation.

No model evaluation or job submission. Default is read-only. Complete archived files and
an inventory remain available; no old result is relabelled with invented provenance.
"""
import argparse
import fcntl
import json
import shutil
from pathlib import Path
import sys
import time
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from redox.core.common import ROOT,UMA,LIBRARY,read_manifest
from redox.qm.uma import model_identity,cache_errors,file_sha256,DEFAULT_MODEL


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--only',required=True,help='comma-separated registered molecule IDs')
    ap.add_argument('--archive-root',type=Path,required=True)
    ap.add_argument('--execute',action='store_true',help='archive stale state directories after acquiring motif locks')
    a=ap.parse_args();ids=set(a.only.split(','));rows=[r for r in read_manifest() if r['id'] in ids]
    if ids-{r['id'] for r in rows}:ap.error('unknown registered ID')
    archive=a.archive_root.resolve()
    if archive==UMA.resolve() or UMA.resolve() in archive.parents:
        ap.error('archives must be outside the canonical UMA tree')
    identity=model_identity(DEFAULT_MODEL);handles=[];moves=[]
    try:
        if a.execute:
            for gid in sorted(ids):
                p=ROOT/'logs/pipeline_locks'/f'{gid}.lock';p.parent.mkdir(parents=True,exist_ok=True)
                h=p.open('a+');handles.append(h)
                fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for row in rows:
            sd=UMA/row['id']/row['state'];rj=sd/'result.json'
            if not sd.exists():continue
            errors=['missing_result'] if not rj.exists() else cache_errors(json.loads(rj.read_text()),row,sd/'relaxed.xyz',
                  sorted((LIBRARY/row['id']/'conformers').glob('conf_*.xyz')),identity,.02)
            if not errors:continue
            files={str(p.relative_to(sd)):file_sha256(p) for p in sd.rglob('*') if p.is_file()}
            target=archive/row['id']/row['state']
            moves.append(dict(id=row['id'],state=row['state'],source=str(sd),archive=str(target),reason=errors,sha256=files))
        if a.execute and moves:
            if archive.exists():ap.error('stale outputs need a fresh archive root; prior archives are preserved')
            archive.mkdir(parents=True)
            inventory=archive/'inventory.json'
            inventory.write_text(json.dumps(dict(status='archive_planned',created=time.time(),states=moves),indent=2)+'\n')
            for m in moves:
                target=Path(m['archive']);target.parent.mkdir(parents=True,exist_ok=True)
                # Move preserves all original bytes; no result field is rewritten.
                shutil.move(m['source'],target)
                found={str(p.relative_to(target)):file_sha256(p) for p in target.rglob('*') if p.is_file()}
                if found!=m['sha256']:raise RuntimeError('archive integrity verification failed')
            inventory.write_text(json.dumps(dict(status='archived_for_recalculation',created=time.time(),states=moves),indent=2)+'\n')
        print(json.dumps(dict(executed=a.execute,states=moves,quantum_jobs_submitted=0),indent=2))
    finally:
        for h in handles:h.close()
if __name__=='__main__':main()
