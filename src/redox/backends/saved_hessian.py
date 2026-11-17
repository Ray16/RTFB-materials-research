"""Reuse one complete Cartesian Hessian only at its declared retained state.

The caller must admit a closed, source-bound producer. The returned matrix is
still unqualified until fresh gradients, sampled derivatives and modes pass.
"""
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path
from types import MethodType

import numpy as np


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@contextmanager
def builder(module,plan,case):
    from redox.core.io import atomic_write_json
    case=Path(case);saved=plan['saved_full_Hessian'];original=module._build_smd_mf;calls=[]
    def build(*args,**kwargs):
        mf=original(*args,**kwargs);factory=mf.Hessian;had_instance='Hessian' in mf.__dict__;previous_instance=mf.__dict__.get('Hessian')
        def hessian(*hargs,**hkwargs):
            # undo_solvent() makes a shallow view of the mean field. An
            # instance override must not propagate into that vacuum view.
            if had_instance:mf.Hessian=previous_instance
            else:delattr(mf,'Hessian')
            try:obj=factory(*hargs,**hkwargs)
            finally:mf.Hessian=hessian
            def kernel(self,*kargs,**kkwargs):
                if kargs or kkwargs or calls:raise ValueError('One full fixed-geometry Hessian checkpoint call required')
                for path,expected in saved['source_sha256'].items():
                    if digest(path)!=expected:raise ValueError('Saved-Hessian backing changed')
                center=json.loads((case/'points/0000/evaluation.json').read_text())
                old=json.loads(Path(saved['center']).read_text());coordinates=np.asarray(center['coordinates_Bohr'])
                gradient=np.asarray(center['gradient_Ha_Bohr']);energy_difference=abs(center['E_Ha']-old['E_Ha'])
                force_difference=float(abs(gradient-np.asarray(old['gradient_Ha_Bohr'])).max())
                if (not center['SCF_converged'] or not math.isfinite(energy_difference) or energy_difference>1e-8
                        or gradient.shape!=(len(plan['symbols']),3) or not np.isfinite(gradient).all() or not math.isfinite(force_difference) or force_difference>1e-6
                        or not np.allclose(coordinates,old['coordinates_Bohr'],rtol=0,atol=1e-12)
                        or not np.allclose(mf.mol.atom_coords(),coordinates,rtol=0,atol=1e-12)
                        or mf.mol.charge!=plan.get('charge',0) or mf.mol.spin!=plan.get('mult',1)-1
                        or mf.xc!='r2scan' or mf.mol.basis!='def2-svp'):
                    raise ValueError('Fresh center does not reproduce the saved matrix state')
                if plan.get('charge',0)!=0 and (not plan.get('charged_saved_Hessian_validation')
                        or saved.get('charge')!=plan['charge'] or saved.get('mult')!=plan.get('mult')
                        or saved.get('source_response_independently_checked') is not True):
                    raise ValueError('Charged matrix reuse requires its exact checked source state')
                matrix=np.load(saved['matrix'],allow_pickle=False);n=len(plan['symbols'])
                if (matrix.dtype!=np.float64 or matrix.shape!=(3*n,3*n) or not np.isfinite(matrix).all()
                        or float(abs(matrix-matrix.T).max())>1e-6):raise ValueError('Complete finite symmetric FP64 source Hessian required')
                calls.append(dict(source_matrix=saved['matrix'],source_matrix_sha256=digest(saved['matrix']),
                    center_energy_difference_Ha=energy_difference,center_gradient_difference_Ha_Bohr=force_difference,
                    full_response_recomputed=False,fresh_derivative_and_minimum_checks_required=True,native_G_qualified=False))
                atomic_write_json(case/'saved_matrix_reuse.json',dict(calls=calls,source=saved,native_G_qualified=False))
                return matrix.reshape(n,3,n,3).transpose(0,2,1,3)
            obj.kernel=MethodType(kernel,obj);return obj
        mf.Hessian=hessian;return mf
    module._build_smd_mf=build
    try:yield calls
    finally:
        module._build_smd_mf=original
        atomic_write_json(case/'saved_matrix_reuse.json',dict(calls=calls,source=saved,native_G_qualified=False))
