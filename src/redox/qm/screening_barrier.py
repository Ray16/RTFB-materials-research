"""Run-local upstream grid-screening barrier correction; no installed-file edits.

Qualification is separate from construction. Changed screening implementations
must not silently reuse or replace source-frozen thermochemistry.
"""
from contextlib import contextmanager
import hashlib
import inspect
from pathlib import Path
import re
import textwrap


def reference_mask(coords,atm,bas,env,block_size,log_cutoff,constants):
    """Independent FP64 per-shell/per-grid-block bound with the same cutoff."""
    import numpy as np
    coords=np.asarray(coords,dtype=np.float64);output=np.zeros(((len(coords)+block_size-1)//block_size,len(bas)),dtype=np.int8)
    for shell,basis in enumerate(bas):
        atom=int(basis[constants['ATOM_OF']]);pointer=int(atm[atom,constants['PTR_COORD']]);center=env[pointer:pointer+3]
        nprim=int(basis[constants['NPRIM_OF']]);exp0=int(basis[constants['PTR_EXP']]);coef0=int(basis[constants['PTR_COEFF']]);angular=int(basis[constants['ANG_OF']])
        for block,start in enumerate(range(0,len(coords),block_size)):
            difference=coords[start:start+block_size]-center;rr=np.sum(difference*difference,axis=1)+1e-300
            bound=np.full(len(rr),1e-300)
            for primitive in range(nprim):bound+=env[coef0+primitive]*np.exp(-env[exp0+primitive]*rr)
            with np.errstate(divide='ignore'):significant=np.log(np.abs(bound))+angular*np.log(rr)/2>log_cutoff
            output[block,shell]=int(significant.any())
    return output


def standalone_source(source,constants,*,corrected=True):
    """Retain the original numerical expressions; add the upstream loop barrier."""
    declaration='__global__\nstatic void _screen_index('
    loop='    for (int grid0 = grid_start; grid0 < grid_stop; grid0 += NG_PER_BLOCK) {\n'
    if source.count(declaration)!=1 or source.count(loop)!=1 or source.count('__syncthreads();')!=1:
        raise ValueError('Exact original screening-kernel source required')
    required={'ATM_SLOTS','PTR_COORD','BAS_SLOTS','ATOM_OF','ANG_OF','NPRIM_OF','PTR_EXP','PTR_COEFF'}
    if set(constants)!=required or any(type(v) is not int or v<0 for v in constants.values()):raise ValueError('Explicit PySCF ABI constants required')
    result=source.replace(declaration,'extern "C" __global__ void redox_screen_index(',1)
    if corrected:result=result.replace(loop,loop+'        __syncthreads(); // upstream PR925: protect previous shared-coordinate readers\n',1)
    definitions='\n'.join('#define '+key+' '+str(value) for key,value in sorted(constants.items()))
    # NVRTC provides CUDA scalar math intrinsics without host C-library
    # headers. Use the exact signed-byte ABI instead of a system stdint.h.
    return 'typedef signed char int8_t;\n#define NG_PER_BLOCK 256\n'+definitions+'\n'+result


def launch(mask,log_cutoff,coords,atm,bas,env,block_size,kernel,stream):
    import numpy as np
    if (mask.dtype!=np.int8 or coords.dtype!=np.float64 or not coords.flags.f_contiguous
            or atm.dtype!=np.int32 or bas.dtype!=np.int32 or env.dtype!=np.float64
            or block_size<256 or block_size%256 or mask.shape!=((len(coords)+block_size-1)//block_size,len(bas))):
        raise ValueError('Screening buffers differ from the installed FP64/int32 ABI')
    kernel((mask.shape[0],(len(bas)+255)//256),(256,),
        (mask,np.float64(log_cutoff),coords,np.int32(len(coords)),np.int32(block_size),atm,np.int32(len(atm)),bas,np.int32(len(bas)),env),stream=stream)


@contextmanager
def corrected_screening(source_path,output):
    import cupy as cp
    from pyscf.gto import mole
    from gpu4pyscf.dft import gen_grid
    from redox.core.io import atomic_write_json
    output=Path(output);output.mkdir(exist_ok=False)
    names=('ATM_SLOTS','PTR_COORD','BAS_SLOTS','ATOM_OF','ANG_OF','NPRIM_OF','PTR_EXP','PTR_COEFF')
    constants={name:int(getattr(mole,name)) for name in names};cuda=standalone_source(Path(source_path).read_text(),constants)
    kernel=cp.RawKernel(cuda,'redox_screen_index',options=('--std=c++11',))
    original=gen_grid.Grids._build_non0ao_idx_cache;body=textwrap.dedent(inspect.getsource(original))
    pattern=r'    libgdft\.GDFTscreen_index\(\n.*?ctypes\.cast\(_env_gpu\.data\.ptr, ctypes\.c_void_p\)\)'
    replacement='    _redox_screening_launch(non0shl_mask, np.log(cutoff), coords, _atm_gpu, _bas_gpu, _env_gpu, block_size, _redox_kernel, stream)'
    changed,count=re.subn(pattern,replacement,body,count=1,flags=re.S)
    if count!=1:raise ValueError('Installed screening host-wrapper template differs')
    (output/'screening_kernel.cu').write_text(cuda);(output/'original_host.py').write_text(body);(output/'corrected_host.py').write_text(changed)
    scope=dict(original.__globals__,_redox_screening_launch=launch,_redox_kernel=kernel);local={}
    exec(compile(changed,str(output/'corrected_host.py'),'exec'),scope,local)
    stats=dict(status='run_local_upstream_screening_correction_active_unqualified',upstream_fix='https://github.com/pyscf/gpu4pyscf/pull/925',constants=constants,
        source_sha256=hashlib.sha256(Path(source_path).read_bytes()).hexdigest(),kernel_sha256=hashlib.sha256(cuda.encode()).hexdigest(),installed_gen_grid_sha256=hashlib.sha256(Path(gen_grid.__file__).read_bytes()).hexdigest(),installed_library_files_modified=False,physics_model_changed=False,implementation_qualified=False,native_G_qualified=False)
    atomic_write_json(output/'implementation.json',stats);gen_grid.Grids._build_non0ao_idx_cache=local['_build_non0ao_idx_cache']
    try:yield stats
    finally:
        gen_grid.Grids._build_non0ao_idx_cache=original
        stats['status']='run_local_screening_context_closed_qualification_separate';atomic_write_json(output/'implementation.json',stats)
