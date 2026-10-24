"""Use the current UMA/DFT entry points with the existing reserved shard launcher."""
from pathlib import Path

from redox.core.common import ROOT
from redox.qm.uma_model_identity import digest
from redox.workflows import launch as original

_commands = original.commands


def commands(module,*args,**kwargs):
    jobs=_commands(module,*args,**kwargs)
    for job in jobs:
        offset=job.index('--')+1
        index=job.index('redox.qm.'+module,offset)
        job[index]='redox.qm.uma_current'
        job[index+1:index+1]=['--stage',module]
    return jobs


def main(argv=None):
    old_commands,writer=original.commands,original.atomic_write_json
    def write(path,record):
        if Path(path).name=='launch.json':
            record['current_entrypoint_source_sha256']={str(p):digest(p) for p in
                (Path(__file__).resolve(),ROOT/'src/redox/qm/uma_current.py',ROOT/'config/uma_execution.json')}
            record['preferred_UMA_backend']='nitro'
        return writer(path,record)
    original.commands,original.atomic_write_json=commands,write
    try:return original.main(argv)
    finally:original.commands,original.atomic_write_json=old_commands,writer


if __name__=='__main__':raise SystemExit(main())
