"""Bounded retries of transient HDF5 checkpoint locks, without bypassing locks."""
from contextlib import contextmanager
import errno
import json
from pathlib import Path
import time

DELAYS=(.1,.2,.5,1.,2.,4.)


def retry_checkpoint(operation,args,kwargs,log,*,sleep=time.sleep,delays=DELAYS):
    for attempt in range(len(delays)+1):
        try:return operation(*args,**kwargs)
        except BlockingIOError as error:
            if error.errno not in (errno.EAGAIN,errno.EWOULDBLOCK) or 'lock' not in str(error).lower():raise
            row=dict(checked=time.time(),operation=operation.__name__,attempt=attempt+1,
                     checkpoint=str(args[1]) if len(args)>1 else None,error=repr(error),
                     action='retry_same_checkpoint_with_locking_enabled' if attempt<len(delays) else 'exhausted_preserved')
            log.parent.mkdir(parents=True,exist_ok=True)
            with log.open('a') as stream:stream.write(json.dumps(row)+'\n')
            if attempt==len(delays):raise
            sleep(delays[attempt])


@contextmanager
def checkpoint_lock_retries(log):
    """Bind only the calling process's PySCF checkpoint functions."""
    from pyscf.scf import chkfile
    original={name:getattr(chkfile,name) for name in ('save_mol','dump_scf')}
    def wrapped(function):
        def call(*args,**kwargs):return retry_checkpoint(function,args,kwargs,Path(log))
        return call
    for name,function in original.items():setattr(chkfile,name,wrapped(function))
    try:yield
    finally:
        for name,function in original.items():setattr(chkfile,name,function)
