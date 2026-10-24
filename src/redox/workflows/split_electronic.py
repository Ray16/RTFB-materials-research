"""Join a checked geometry and a canonical energy with explicit backend identity."""
import copy
import math

from redox.core.protocol import ACTIVE_SP


def composite_result(geometry,active,geometry_sha1,energy_record_path,*,backend):
    from redox.qm.sp import _ENERGY_KEYS,_adoptable
    if backend not in ('cpu','gpu'):
        raise ValueError('An explicit supported electronic backend is required')
    if any(geometry.get(k) is not True for k in ('optimized','geometry_optimization_converged','optimization_SCF_converged')):
        raise ValueError('A converged optimizer/SCF parent is required')
    if (active.get('source')!='computed:'+backend or active.get('protocol')!=ACTIVE_SP
            or active.get('geom_sha1')!=geometry_sha1
            or any(type(active.get(k)) is not int or type(geometry.get(k)) is not int for k in ('charge','mult'))
            or any(active.get(k)!=geometry.get(k) for k in ('charge','mult'))):
        raise ValueError('The active energy differs in backend, protocol, state or geometry')
    if (any(active.get(k) is not True for k in ('converged_smd','converged_gas'))
            or any(type(active.get(k)) not in (int,float) or not math.isfinite(active[k]) for k in ('e_smd_eV','e_gas_eV'))):
        raise ValueError('Both finite, converged active energies are required')
    result=copy.deepcopy(geometry)
    result.update({key:active.get(key) for key in _ENERGY_KEYS})
    result.update(geometry_only_result=False,split_execution=True,electronic_backend=backend,
        geometry_backend=geometry['backend'],electronic_record_source=str(energy_record_path),
        energy_geom_sha1=geometry_sha1,energy_protocol=copy.deepcopy(ACTIVE_SP),
        sp_xc=ACTIVE_SP['xc'],sp_basis=ACTIVE_SP['basis'],sp_nlc=ACTIVE_SP['nlc'],
        sp_disp=ACTIVE_SP['disp'],solvent=ACTIVE_SP['solvent'])
    if not _adoptable(result,charge=geometry['charge'],mult=geometry['mult'],geom_hash=geometry_sha1):
        raise ValueError('The exact energy cannot be adopted by the normal native stage')
    return result
