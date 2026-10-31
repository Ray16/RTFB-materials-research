"""Read-only orbital checkpoints and transparent response-call instrumentation.

These helpers do not supply thermochemistry or qualify a restarted Hessian.
"""
from contextlib import contextmanager
import faulthandler
import hashlib
import io
import json
import math
import os
from pathlib import Path
import threading
import time


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_array(path,value):
    """Publish one durable checkpoint after its complete NPY write."""
    import numpy as np
    path=Path(path)
    if path.exists():raise FileExistsError('Preserve earlier response arrays')
    array=np.asarray(value.get() if hasattr(value,'get') else value)
    if array.dtype!=np.float64 or not np.isfinite(array).all():raise ValueError('Finite FP64 response checkpoint required')
    temporary=path.with_name(path.name+'.partial')
    with temporary.open('xb') as stream:np.save(stream,array,allow_pickle=False);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)
    return dict(path=str(path),shape=list(array.shape),dtype='float64',file_sha256=digest(path),array_sha256=hashlib.sha256(array.tobytes()).hexdigest(),native_G_qualified=False)


def contiguous_payload(handle,name,content):
    """Read only native FP64 contiguous payloads with verified in-file bounds.

    The installed HDF5 high-level read stalls on the real source checkpoint.
    Metadata offsets avoid that path; orbital SHA-256 checks below remain exact.
    Chunked, compressed, external or unallocated datasets are unsupported.
    """
    import h5py
    import numpy as np
    dataset=handle[name];creation=dataset.id.get_create_plist();offset=dataset.id.get_offset()
    if (dataset.dtype!=np.dtype('float64') or creation.get_layout()!=h5py.h5d.CONTIGUOUS
            or creation.get_nfilters()!=0 or creation.get_external_count()!=0
            or offset is None or offset<0):raise ValueError('Only allocated native FP64 contiguous checkpoint datasets are supported')
    count=math.prod(dataset.shape);size=count*8
    if offset+size>len(content):raise ValueError('Contiguous dataset extends beyond the source checkpoint')
    return np.frombuffer(content,dtype=np.float64,count=count,offset=offset).reshape(dataset.shape)


def read_checkpoint(raw):
    import h5py
    import numpy as np
    raw=Path(raw);identity_path=raw/'execution/rks_make_h1_state_identity.json'
    evaluation_path=raw/'SCF_evaluations/0000/evaluation.json';settings_path=evaluation_path.with_name('settings.json')
    checkpoint=evaluation_path.with_name('SCF.h5')
    captured={path:path.read_bytes() for path in (identity_path,evaluation_path,settings_path,checkpoint,raw/'geometry.xyz')}
    identity=json.loads(captured[identity_path]);evaluation=json.loads(captured[evaluation_path]);settings=json.loads(captured[settings_path])
    matrix=raw/'execution/rks_make_h1_0.npy'
    expected=dict(charge=0,mult=1,xc='r2scan',SCF_converged=True,basis="'def2-svp'",stage='rks_make_h1')
    if any(identity.get(k)!=v for k,v in expected.items()) or evaluation.get('SCF_converged') is not True:
        raise ValueError('A converged exact neutral r2SCAN/SVP orbital checkpoint is required')
    if any(settings.get(k)!=v for k,v in dict(charge=0,mult=1,basis='def2-svp',xc='r2scan',disp='d4',solvent='acetonitrile',grid_level=3).items()):
        raise ValueError('The source checkpoint numerical/solvent settings differ')
    if np.max(abs(np.asarray(identity['coordinates_Bohr'])-np.asarray(evaluation['coordinates_Bohr'])))>1e-12:
        raise ValueError('Orbital and source evaluation geometries differ')
    matches=[v for p,v in identity['arrays'].items() if Path(p).name==matrix.name]
    if len(matches)!=1 or digest(matrix)!=matches[0]:raise ValueError('First-derivative matrix is not bound to the orbital state')
    arrays={}
    with h5py.File(io.BytesIO(captured[checkpoint]),'r') as handle:
        energy=float(contiguous_payload(handle,'scf/e_tot',captured[checkpoint]))
        for name in ('mo_coeff','mo_energy','mo_occ'):
            array=contiguous_payload(handle,'scf/'+name,captured[checkpoint]);record=identity[name]
            if (array.dtype!=np.float64 or list(array.shape)!=record['shape'] or record['dtype']!='float64'
                    or not np.isfinite(array).all() or hashlib.sha256(array.tobytes()).hexdigest()!=record['sha256']):
                raise ValueError('Saved SCF orbitals do not exactly match the derivative frame')
            arrays[name]=array
    if not np.isfinite(energy) or abs(energy-identity['SCF_E_Ha'])>1e-12 or abs(energy-evaluation['E_Ha'])>1e-12:
        raise ValueError('SCF checkpoint energy differs from its source evaluation')
    h1=np.load(matrix,mmap_mode='r',allow_pickle=False);n=len(identity['coordinates_Bohr'])
    if h1.dtype!=np.float64 or h1.shape!=(n,3,arrays['mo_coeff'].shape[1],int((arrays['mo_occ']>0).sum())) or not np.isfinite(h1).all():
        raise ValueError('The saved first-derivative tensor has invalid shape or values')
    hashes={str(path):hashlib.sha256(content).hexdigest() for path,content in captured.items()};hashes[str(matrix)]=matches[0]
    if any(digest(path)!=value for path,value in hashes.items()):raise ValueError('A source checkpoint changed during inspection')
    return dict(identity=identity,settings=settings,energy_Ha=energy,source_sha256=hashes,
        reader='hash-checked contiguous native-FP64 payload; HDF5 metadata only'),arrays,h1


@contextmanager
def response_memory_cap(module,cap_bytes):
    """Reduce automatic RHS batching; never advertise more than actual memory."""
    if cap_bytes!=3*2**30:raise ValueError('This diagnostic profile uses a 3-GiB available-memory ceiling')
    original=module.get_avail_mem
    module.get_avail_mem=lambda *args,**kwargs:min(original(*args,**kwargs),cap_bytes)
    try:yield
    finally:module.get_avail_mem=original


@contextmanager
def trace_solver(module,directory,*,trace_seconds=30,checkpoint_root=None,checkpoint_identity=None):
    """Keep the frozen solver and response values; add host-side event timing."""
    from redox.core.io import atomic_write_json
    directory=Path(directory);directory.mkdir(exist_ok=False)
    if checkpoint_root is not None:Path(checkpoint_root).mkdir(exist_ok=False)
    original=module.residual_krylov;events=[];calls=[]
    def event(kind,**fields):
        row=dict(kind=kind,monotonic=time.monotonic(),UTC=time.time(),**fields)
        events.append(row)
        with (directory/'events.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
        atomic_write_json(directory/'latest.json',row)
    def solve(aop,b,*args,**kwargs):
        import numpy as np
        index=len(calls);calls.append(dict(index=index,shape=list(b.shape),status='entered'))
        rhs_sha=None
        if checkpoint_root is not None:
            rhs=np.asarray(b.get() if hasattr(b,'get') else b);rhs_sha=hashlib.sha256(rhs.tobytes()).hexdigest();del rhs
        event('solver_entered',index=index,shape=list(b.shape));operator_calls=0
        def operator(vector):
            nonlocal operator_calls
            iteration=operator_calls;operator_calls+=1;started=time.monotonic()
            event('operator_entered',solver=index,operator=iteration,shape=list(vector.shape))
            value=aop(vector)
            # GPU work can be asynchronous. This is host-call wall time only;
            # instrumentation introduces no synchronization or array arithmetic.
            event('operator_returned',solver=index,operator=iteration,host_elapsed_seconds=time.monotonic()-started)
            return value
        result=original(operator,b,*args,**kwargs)
        if checkpoint_root is not None:
            checkpoint=save_array(Path(checkpoint_root)/f'block_{index:03d}.npy',result)
            atomic_write_json(Path(checkpoint_root)/f'block_{index:03d}.json',dict(checkpoint=checkpoint,rhs_sha256=rhs_sha,identity=checkpoint_identity,
                original_solver_source_sha256=digest(original.__code__.co_filename),kind='raw_Krylov_solution_before_occupied_block_reconstruction',native_G_qualified=False))
            calls[index]['checkpoint']=checkpoint
        calls[index].update(status='returned',operator_calls=operator_calls)
        event('solver_returned',index=index,operator_calls=operator_calls)
        return result
    with (directory/'stacks.log').open('x') as stacks:
        # The C watchdog used by dump_traceback_later walks changing Python
        # thread/frame state without the GIL (CPython issue 116008). A normal
        # Python thread enters dump_traceback with the GIL held instead.
        stopped=threading.Event()
        def sample():
            while not stopped.wait(trace_seconds):
                faulthandler.dump_traceback(file=stacks,all_threads=True)
        sampler=threading.Thread(target=sample,name='response-GIL-held-stack-sampler',daemon=True);sampler.start()
        module.residual_krylov=solve
        try:yield calls
        finally:
            module.residual_krylov=original;stopped.set();sampler.join()
            atomic_write_json(directory/'trace.json',dict(calls=calls,event_count=len(events),solver_replaced=False,
                array_values_modified_by_instrumentation=False,GPU_synchronization_added=checkpoint_root is not None,asynchronous_C_watchdog_used=False,native_G_qualified=False))
