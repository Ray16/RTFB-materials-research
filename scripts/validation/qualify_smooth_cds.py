#!/usr/bin/env python
"""CPU component qualification of an explicitly experimental smooth CDS model.

The installed PySCF experimental surface model is compared to the unchanged
legacy CDS implementation. This is a method-sensitivity study, never an energy,
geometry or free-energy substitution in the primary candidate records.
"""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def gradient_functions(energy_module, snapshot):
    """Use the installed derivative formulas with their intended energy tables.

    The experimental gradient's imports reference names absent from the legacy
    module. Load only its function definitions into a private namespace with
    the matching experimental tables. Installed files are never modified.
    """
    from pyscf.data import radii
    from pyscf.solvent import pcm
    from pyscf.solvent.grad import pcm as pcm_grad
    from pyscf.lib import logger

    source = Path(energy_module.__file__).parent / 'grad/smd_experiment.py'
    tree = ast.parse(source.read_text())
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = dict(vars(energy_module), np=np, radii=radii, pcm=pcm,
                     pcm_grad=pcm_grad, smd=energy_module, logger=logger)
    snapshot.write_bytes(source.read_bytes())
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), 'exec'), namespace)
    return namespace['get_cds'], source


def solvent_object(xyz, charge=0, mult=1, grid=590):
    from ase.io import read
    from pyscf import gto

    atoms = read(xyz)
    if set(atoms.get_chemical_symbols()) - {'C', 'H', 'N', 'O', 'Cl'}:
        raise ValueError('This pilot is scoped to C/H/N/O/Cl states; no metal extrapolation')
    mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions)),
                unit='Angstrom', charge=charge, spin=mult-1,
                basis='def2-svp', verbose=0)
    return SimpleNamespace(mol=mol, solvent='acetonitrile', _solvent='acetonitrile',
                           solvent_descriptors=None, sasa_ng=grid)


def sampled_gradient_check(obj, energy, gradient, selectors):
    original = obj.mol.atom_coords().copy()
    analytic = np.asarray(gradient(obj))
    if analytic.shape != original.shape or not np.isfinite(analytic).all():
        raise ValueError('Invalid smooth-CDS gradient')
    rows = []
    try:
        for atom, axis in selectors:
            for step in (1e-4, 5e-4, 1e-3):
                plus, minus = original.copy(), original.copy()
                plus[atom, axis] += step
                minus[atom, axis] -= step
                obj.mol.set_geom_(plus, unit='Bohr')
                ep = energy(obj)
                obj.mol.set_geom_(minus, unit='Bohr')
                em = energy(obj)
                numerical = (ep-em)/(2*step)
                rows.append(dict(atom=atom, axis=axis, step_Bohr=step,
                    analytic_Ha_Bohr=float(analytic[atom, axis]),
                    numerical_Ha_Bohr=float(numerical),
                    absolute_error_Ha_Bohr=float(abs(numerical-analytic[atom, axis]))))
    finally:
        obj.mol.set_geom_(original, unit='Bohr')
    return dict(tested_components=len(selectors), total_components=original.size,
                checks=rows, maximum_error_Ha_Bohr=max(r['absolute_error_Ha_Bohr'] for r in rows))


def contact_scan(obj, energy):
    original = obj.mol.atom_coords().copy()
    bohr = .529177210903
    direction = original[34]-original[19]
    direction /= np.linalg.norm(direction)
    values = {}
    try:
        for offset in (-1e-3, -1e-4, -4e-5, -2e-5, -1e-5, 0.,
                       1e-5, 2e-5, 4e-5, 1e-4, 1e-3):
            xyz = original.copy()
            xyz[34] = xyz[19]+direction*((4.2+offset)/bohr)
            obj.mol.set_geom_(xyz, unit='Bohr')
            values[offset] = float(energy(obj))
    finally:
        obj.mol.set_geom_(original, unit='Bohr')
    def slope(side):
        x = np.array([side*v/bohr for v in (1e-5, 2e-5, 4e-5)])
        y = np.array([values[side*v] for v in (1e-5, 2e-5, 4e-5)])
        return float(np.dot(x-x.mean(), y-y.mean())/np.dot(x-x.mean(), x-x.mean()))
    left, right = slope(-1), slope(1)
    return dict(one_sided_slope_jump_Ha_Bohr=right-left,
        centered_curvature_Ha_Bohr2={str(h):(values[h]-2*values[0.]+values[-h])/(h/bohr)**2
                                   for h in (1e-3, 1e-4, 1e-5)},
        energy_Ha_by_offset_A={str(k):v for k,v in values.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--contact', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    from redox.workflows import polymer_campaign as api
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=8)
        os.sched_setaffinity(0, cpu_allocation(1))
    from pyscf.solvent import smd as legacy
    from pyscf.solvent import smd_experiment as smooth

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    grad, grad_source = gradient_functions(smooth, out/'experimental_gradient_source.py')
    report = json.loads(args.source.read_text())
    api.verify_sources(report['source_sha256'])
    contact = json.loads(args.contact.read_text())
    api.verify_sources(contact['source_sha256'])
    states = {'contact': dict(state_dir=str(Path(contact['source_geometry']).parent))}
    states.update({k:dict(v) for k,v in report['reagent_states'].items()})
    states['empty'] = report['empty_reference']
    states['capped'] = report['local_capped_reference']
    states.update({'site'+str(r['site']):r for r in report['singles']})
    sources = {str(p.resolve()):digest(p) for p in (Path(__file__), args.source, args.contact,
        Path(legacy.__file__), Path(smooth.__file__), grad_source)}
    rows = {}
    record = dict(status='experimental_smooth_CDS_component_study_running', started=time.time(),
        CPU_affinity=sorted(os.sched_getaffinity(0)), states=rows,
        acceptance=dict(sampled_gradient_maximum_Ha_Bohr=1e-6,
            finest_grid_energy_change_kJ_mol=.5, rigid_rotation_spread_kJ_mol=.5),
        method_changed=True, primary_records_changed=False, GPU_jobs_started=0,
        source_sha256=sources, full_DFT_or_G_qualified=False)
    atomic_write_json(out/'summary.json',record)
    try:
        for name, state in states.items():
            directory = Path(state['state_dir'])
            xyz = directory/'opt.xyz'
            raw = json.loads((directory/'result.json').read_text())
            obj = solvent_object(xyz, raw['charge'], raw['mult'])
            sources[str(xyz)] = digest(xyz)
            points = []
            for grid in (590, 974, 2030):
                obj.sasa_ng = grid
                points.append(dict(grid=grid, smooth_E_Ha=float(smooth.get_cds(obj))))
            legacy_energy = float(legacy.get_cds_legacy(obj)[0])
            values = dict(grid_points=points, legacy_E_Ha=legacy_energy,
                fixed_geometry_model_shift_kJ_mol=(points[-1]['smooth_E_Ha']-legacy_energy)*2625.4996394799,
                finest_grid_change_kJ_mol=abs(points[-1]['smooth_E_Ha']-points[-2]['smooth_E_Ha'])*2625.4996394799)
            if name == 'contact':
                values['contact_scan'] = contact_scan(obj, smooth.get_cds)
                values['sampled_gradient'] = sampled_gradient_check(obj, smooth.get_cds, grad,
                    [(19,0),(34,0),(0,1),(obj.mol.natm-1,2)])
                original = obj.mol.atom_coords().copy()
                energies = [float(smooth.get_cds(obj))]
                rng = np.random.default_rng(20261025)
                for _ in range(3):
                    rotation, _ = np.linalg.qr(rng.normal(size=(3,3)))
                    obj.mol.set_geom_((original-original.mean(axis=0))@rotation, unit='Bohr')
                    energies.append(float(smooth.get_cds(obj)))
                obj.mol.set_geom_(original, unit='Bohr')
                values['rigid_rotation_spread_kJ_mol'] = (max(energies)-min(energies))*2625.4996394799
                # Establish equality of the CPU legacy reference and the exact
                # prior contact scan before interpreting model differences.
                legacy_scan = contact_scan(obj, lambda o:legacy.get_cds_legacy(o)[0])
                prior = next(r['E_CDS_Ha'] for r in contact['points'] if r['offset_A']==0.)
                error = abs(legacy_scan['energy_Ha_by_offset_A']['0.0']-prior)
                values['CPU_legacy_contact_reference_error_Ha'] = error
                if error>1e-11:
                    raise ValueError('CPU legacy reference differs from prior exact contact scan')
            rows[name] = values
            atomic_write_json(out/'summary.json',record)
        atom = rows['contact']
        record['component_checks_passed'] = bool(
            atom['sampled_gradient']['maximum_error_Ha_Bohr']<=1e-6
            and atom['rigid_rotation_spread_kJ_mol']<=.5
            and max(v['finest_grid_change_kJ_mol'] for v in rows.values())<=.5)
        record['status'] = 'experimental_smooth_CDS_component_study_complete'
        record['interpretation'] = ('Experimental model sensitivity only. Passing component checks does not '
            'qualify a full solvent protocol, geometry, native Hessian or free energy.')
        api.verify_sources(sources)
    except BaseException as error:
        record.update(status='component_study_failed_preserved',error=repr(error))
        raise
    finally:
        record['finished'] = time.time()
        atomic_write_json(out/'summary.json',record)
    print(json.dumps(dict(status=record['status'],component_checks_passed=record['component_checks_passed'])))


if __name__ == '__main__':
    main()
