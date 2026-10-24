"""Join separately qualified geometry and active energies without losing parents."""
import copy
import math

from redox.core.protocol import ACTIVE_SP


def composite_result(geometry, active, geometry_sha1, energy_record_path):
    """Materialize a new state's two-stage record; never change a parent record.

    The canonical SP record remains the authority for energy. The explicit
    combined result lets the native stage adopt that same energy instead of
    repeating a large hybrid calculation on its Hessian GPU.
    """
    from redox.qm.sp import _ENERGY_KEYS, _adoptable

    if (geometry.get('geometry_optimization_converged') is not True
            or geometry.get('optimization_SCF_converged') is not True
            or geometry.get('optimized') is not True):
        raise ValueError('A converged geometry parent is required')
    if (active.get('protocol') != ACTIVE_SP or active.get('geom_sha1') != geometry_sha1
            or active.get('charge') != geometry.get('charge')
            or active.get('mult') != geometry.get('mult')
            or active.get('source') != 'computed:cpu'):
        raise ValueError('The CPU active energy differs in geometry, state or protocol')
    if any(active.get(k) is not True for k in ('converged_smd', 'converged_gas')):
        raise ValueError('Both gas and SMD SCF calculations must converge')
    if any(type(active.get(k)) not in (int, float) or not math.isfinite(active[k])
           for k in ('e_smd_eV', 'e_gas_eV')):
        raise ValueError('Finite active energies required')
    result = copy.deepcopy(geometry)
    result.update({key: active.get(key) for key in _ENERGY_KEYS})
    result.update(geometry_only_result=False, split_execution=True,
        electronic_backend='cpu', geometry_backend=geometry['backend'],
        electronic_record_source=str(energy_record_path),
        energy_geom_sha1=geometry_sha1, energy_protocol=copy.deepcopy(ACTIVE_SP),
        sp_xc=ACTIVE_SP['xc'], sp_basis=ACTIVE_SP['basis'],
        sp_nlc=ACTIVE_SP['nlc'], sp_disp=ACTIVE_SP['disp'], solvent=ACTIVE_SP['solvent'])
    if not _adoptable(result, charge=geometry['charge'], mult=geometry['mult'], geom_hash=geometry_sha1):
        raise ValueError('The native stage cannot reuse the exact CPU parent normally')
    return result
