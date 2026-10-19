"""Private, strict native-level relaxation of two independently prepared basins."""
import argparse
import json
import os
from pathlib import Path
import socket
import time
import traceback

from redox.core.common import read_manifest
from redox.core.io import atomic_write_json
from redox.qm.uma import file_sha256
from redox.workflows.resources import verify_gate_owner


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--task',type=int,required=True)
    args=parser.parse_args();base=Path(__file__).resolve().parent
    plan=json.loads((base/'native_seed_plan.json').read_text());task=plan['tasks'][args.task]
    verify_gate_owner()
    assert len(os.sched_getaffinity(0))==2
    assert all(file_sha256(p)==h for p,h in plan['source_sha256'].items())
    from redox.qm.native_thermal import NATIVE_THERMAL,implementation_controls,ensure_native_thermal,select_policy
    implementation_controls()
    output=base/'native_seeds'/task['key'];output.mkdir(parents=True,exist_ok=False)
    work=output/'work';work.mkdir()
    report=dict(status='computing',stage='native_level_optimization',task=task,
        host=socket.gethostname(),pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),
        canonical_outputs_modified=False,production_admission=False,
        lower_unqualified_basin_not_replaced=True)
    def save():atomic_write_json(output/'status.json',report)
    save()
    try:
        from ase.io import read,write
        from pyscf import gto
        from pyscf.geomopt.geometric_solver import kernel
        from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
        from redox.qm.sp import ensure
        from redox.core.protocol import load_record
        from redox.validation.protonation import hydrogen_topology_status
        from redox.qm.uma import geometry_checks
        atoms=read(task['seed']);q,m=plan['charge'],plan['mult'];symbols=atoms.get_chemical_symbols()
        row=next(r for r in read_manifest() if r['id']==plan['id'] and r['state']==plan['state'])
        basis=NATIVE_THERMAL['anion_basis'];write(work/'seed.xyz',atoms)
        mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),charge=q,spin=m-1,
            basis=basis,ecp=_ecp_for(symbols,basis),verbose=0,max_memory=6000)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
        mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
        mf.kernel();assert mf.converged,'Initial SCF did not converge'
        gradient=mf.nuc_grad_method();gradient.grid_response=True;scanner=gradient.as_scanner()
        converged,finalmol=kernel(scanner,maxsteps=120,assert_convergence=True,
            convergence_energy=1e-8,convergence_grms=3e-5,convergence_gmax=5e-5,
            convergence_drms=1e-4,convergence_dmax=2e-4)
        assert converged and scanner.converged,'Strict geometry/SCF convergence failed'
        atoms.positions=finalmol.atom_coords(unit='Angstrom');write(work/'opt.xyz',atoms)
        geometry_energy_Ha=float(scanner.e_tot)
        del mf,gradient,scanner;_free_gpu_pool()
        assert not geometry_checks(row['smiles'],work/'opt.xyz',q,m),'Geometry integrity failed'
        assert hydrogen_topology_status(row['smiles'],work/'opt.xyz')=='ok','Hydrogen topology changed'
        raw=dict(id=plan['id'],state=plan['state'],charge=q,mult=m,optimized=True,
            geometry_optimization_converged=True,optimization_SCF_converged=True,
            opt_gradient_grid_response=True,opt_xc='r2scan',opt_basis=basis,opt_disp='d4',
            solvent='acetonitrile',opt_out=str(work/'opt.xyz'),g_thermal_eV=None,
            distinct_basin_seed=task['seed'],distinct_basin_seed_sha256=file_sha256(task['seed']),
            canonical_outputs_modified=False)
        atomic_write_json(work/'result.json',raw)
        report.update(stage='active_energy',geometry_level_energy_Ha=geometry_energy_Ha);save()
        _,energy=ensure(plan['id'],plan['state'],work,q,m,backend='gpu')
        assert load_record(work,q,m)
        report.update(stage='native_thermal',active_e_smd_eV=energy['e_smd_eV']);save()
        select_policy(work)
        _,thermal=ensure_native_thermal(plan['id'],plan['state'],state_dir=work,memory_cap_gib=24)
        assert all(file_sha256(p)==h for p,h in plan['source_sha256'].items())
        report.update(status='private_native_state_verified_pending_conformer_and_energy_review',
            thermal_eV=thermal['g_thermal_total_eV'],minimum_frequency_cm=min(thermal['frequencies_cm']),
            output_sha256={str(p):file_sha256(p) for p in work.rglob('*') if p.is_file()})
    except BaseException as exc:
        report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        report['finished']=time.time();save()


if __name__=='__main__':main()
