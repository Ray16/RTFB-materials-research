"""Isolated JSON-lines adapter to the qualified xTB thermodynamic kernel.

Only standard-library modules are imported. This process never imports PySCF,
D4 Python bindings, NumPy or a GPU backend, avoiding incompatible SONAME reuse.
"""
import argparse
import ctypes as C
import hashlib
import json
import math
import os
from pathlib import Path
import sys


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--sha256',required=True)
    args=parser.parse_args()
    actual=hashlib.sha256(args.library.read_bytes()).hexdigest()
    if actual!=args.sha256:raise RuntimeError('Unqualified xTB library')
    library=C.CDLL(str(args.library))
    mass=getattr(library,'__xtb_param_atomicmass_MOD_getatomicmassnumber')
    mass.argtypes=[C.POINTER(C.c_int)];mass.restype=C.c_double
    kernel=getattr(library,'__xtb_thermo_MOD_thermodyn')
    integer,real=C.POINTER(C.c_int),C.POINTER(C.c_double)
    kernel.argtypes=[integer,real,real,real,real,integer,integer,real,real,real,
                     integer,real,real,real,real,real,real,real,real,integer]
    kernel.restype=None
    metadata=dict(pid=os.getpid(),host=os.uname().nodename,
                  process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                  cpu_affinity=sorted(os.sched_getaffinity(0)),interpreter=sys.executable,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  library_sha256=actual,
                  thread_limits={k:os.environ.get(k) for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')})
    for line in sys.stdin:
        request=json.loads(line)
        try:
            if request['action']=='masses':
                numbers=request['numbers']
                if any(type(z)!=int or not 1<=z<=118 for z in numbers):raise ValueError('invalid atomic number')
                result=[mass(C.byref(C.c_int(z))) / (1.660539040e-27/9.10938356e-31) for z in numbers]
            elif request['action']=='thermal':
                p=request['parameters'];frequencies=p['frequencies_cm']
                if any(not math.isfinite(v) or v<=0 for v in frequencies):raise ValueError('invalid internal frequency')
                values=[v/219474.63067 for v in frequencies]
                vib=(C.c_double*len(values))(*values)
                aa,bb,cc=[C.c_double(float(v)) for v in p['rotational_constants_cm']]
                average,mass_total,symmetry,temperature,cutoff=[C.c_double(float(p[k])) for k in
                    ('average_moment_SI','total_mass_amu','symmetry_number','temperature_K','cutoff_cm')]
                unit,linear,atom,count,printing=[C.c_int(v) for v in (6,int(p['linear']),int(p['atom']),len(values),0)]
                escf,energy,enthalpy,free,entropy,zp=[C.c_double(v) for v in (0.,0.,0.,0.,0.,.5*sum(values))]
                ref=C.byref
                kernel(ref(unit),ref(aa),ref(bb),ref(cc),ref(average),ref(linear),ref(atom),
                       ref(symmetry),ref(mass_total),vib,ref(count),ref(escf),ref(temperature),ref(cutoff),
                       ref(energy),ref(enthalpy),ref(free),ref(entropy),ref(zp),ref(printing))
                result=dict(G_Ha=free.value,H_Ha=energy.value,thermal_enthalpy_Ha=enthalpy.value,
                            TS_Ha=entropy.value,ZPE_Ha=zp.value)
                if not all(math.isfinite(v) for v in result.values()):raise RuntimeError('nonfinite kernel output')
            else:raise ValueError('unknown kernel action')
            response=dict(id=request['id'],status='ok',result=result,worker=metadata)
        except Exception as exc:
            response=dict(id=request.get('id'),status='error',error=repr(exc),worker=metadata)
        print(json.dumps(response,allow_nan=False),flush=True)


if __name__=='__main__':main()
