#!/usr/bin/env python
"""Isolated neutral-SVPD optimization/Hessian controls; never primary admission.

Primary records and their certificates are not overwritten or relabeled. The
alternative surface has its own protocol identity and independently tested
Hessian columns. Only the basis-independent, certified projected qRRHO integrator
is reused. Comparisons hold the existing anion and Fc reference fixed.
"""
import argparse,hashlib,json,math,os,time
from pathlib import Path
import numpy as np
from redox.core.common import ROOT,DFT
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner

PANEL=('bm_1_4_benzoquinone','bm_1_2_naphthoquinone','bm_9_10_phenanthrenequinone',
       'bm_n_methylphthalimide','bm_n_n_diisopropylpyromellitic_diimide')
SETTINGS=dict(xc='r2scan',disp='d4',basis='def2-svpd',solvent='acetonitrile',solvent_model='SMD',
    grid_level=3,grid_response=True,conv_tol=1e-12,conv_tol_grad=1e-8,
    force_acceptance_Ha_Bohr=.00045,Hessian_asymmetry_acceptance=1e-6,
    Hessian_column_relative_acceptance=.01,min_frequency_cm=-5,temperature_K=298.15)
OPTIMIZER=dict(convergence_energy=1e-8,convergence_grms=3e-5,convergence_gmax=5e-5,
    convergence_drms=1e-4,convergence_dmax=2e-4)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(campaign,key):
    verify_gate_owner()
    plan=json.loads((campaign/'plan.json').read_text())
    if any(sha(p)!=h for p,h in plan['source_sha256'].items()):raise ValueError('Control sources changed')
    task=next(t for t in plan['tasks'] if t['key']==key)
    if sha(task['input'])!=task['input_sha256']:raise ValueError('Control input changed')
    if len(os.sched_getaffinity(0))!=2 or any(os.environ.get(k)!='2' for k in (
        'OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two disjoint CPU slots and all native caps=2 required')
    out=Path(plan['raw_root'])/key;out.mkdir(parents=True,exist_ok=False)
    status=dict(status='running',task=task,started=time.time(),host=os.uname().nodename,pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)),primary_native_admission=False)
    atomic_write_json(out/'status.json',status)
    try:
        from redox.qm import dft
        from redox.workflows.geometry_checkpoints import checkpoint_callback
        from redox.validation.protonation import hydrogen_topology_status
        from redox.qm.sp import ensure
        from redox.core.protocol import load_record
        from redox.qm.projected_thermal import ensure_projected_thermal
        from redox.qm.hessian_memory import execution_budget,stage_budget_bytes
        from redox.qm.native_thermal import NATIVE_THERMAL,installed_versions
        from ase.io import read
        from pyscf import gto
        # Explicit process-local strict optimizer settings; primary producers unchanged.
        previous=dft._build_smd_mf
        def build(*args,**kwargs):
            mf=previous(*args,**kwargs)
            if args[1]=='r2scan':
                mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
            return mf
        dft._build_smd_mf=build
        result=dft.dft_smd(Path(task['input']),0,1,opt_basis='def2-svpd',opt_basis_anion='def2-svpd',
            backend='gpu',do_freq=False,opt_out=out/'opt.xyz',max_opt_steps=200,
            conv_params=dict(OPTIMIZER,callback=checkpoint_callback(out/'geometry_checkpoints')))
        if not result.get('optimized') or hydrogen_topology_status(task['smiles'],out/'opt.xyz')!='ok':
            raise ValueError('Alternative geometry failed convergence or exact topology')
        result.update(id=key,state='neu',n_e=0,control_settings=SETTINGS)
        atomic_write_json(out/'result.json',result)
        ensure(key,'neu',out,0,1,backend='gpu');energy=load_record(out,0,1)
        if energy is None:raise ValueError('Active SP failed ordinary energy checks')
        atoms=read(out/'opt.xyz');symbols=atoms.get_chemical_symbols();coords=atoms.positions/.52917721067
        memory=execution_budget(NATIVE_THERMAL,6,False)
        import gpu4pyscf.df.hessian.rhf as dr
        import gpu4pyscf.df.hessian.uhf as du
        for module in (dr,du):
            if not hasattr(module,'get_avail_mem'):continue
            old=module.get_avail_mem
            module.get_avail_mem=lambda *a,_old=old,_name=module.__name__,**kw:min(_old(*a,**kw),stage_budget_bytes(memory,_name,__import__('sys')._getframe(1).f_code.co_name))
        def cpu(x):return np.asarray(x.get() if hasattr(x,'get') else x)
        def solve(x):
            mol=gto.M(atom=list(zip(symbols,x.tolist())),unit='Bohr',charge=0,spin=0,basis='def2-svpd',verbose=0,max_memory=6000)
            mf=build(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.kernel()
            if not mf.converged:raise ValueError('Alternative SCF failed')
            return mol,mf
        mol,mf=solve(coords);diag=dft._orbital_diag(mf,0)
        grad=mf.nuc_grad_method();grad.grid_response=True;gradient=cpu(grad.kernel())
        if abs(gradient).max()>.00045 or abs(diag['spin_contam'])>.1:raise ValueError('Alternative force/spin gate failed')
        mf.cphf_grids=mf.grids.copy();hobj=mf.Hessian();hobj.grid_response=True;hobj.max_cycle=100
        H=cpu(hobj.kernel());matrix=H.transpose(0,2,1,3).reshape(3*len(atoms),3*len(atoms))
        native=out/'control_surface.json';work=native.with_suffix('');work.mkdir()
        np.save(work/'hessian_cartesian.npy',matrix)
        asym=float(abs(matrix-matrix.T).max())
        if not np.isfinite(matrix).all() or asym>1e-6:raise ValueError('Alternative Hessian finite/symmetry gate failed')
        del mf,hobj,grad;dft._free_gpu_pool();checks=[]
        oxygen=next((i for i,s in enumerate(symbols) if s=='O'),len(atoms)-1)
        for atom,axis in ((0,1),(oxygen,2)):
            for step in (.0005,.001):
                values=[]
                for sign in (-1,1):
                    x=coords.copy();x[atom,axis]+=sign*step
                    _,fm=solve(x);g=fm.nuc_grad_method();g.grid_response=True;values.append(cpu(g.kernel()))
                    del fm,g;dft._free_gpu_pool()
                numerical=(values[1]-values[0])/(2*step)
                error=float(np.linalg.norm(numerical-H[atom,:,axis,:])/np.linalg.norm(numerical))
                checks.append(dict(atom_index=atom,axis=axis,step_Bohr=step,relative_error=error))
                atomic_write_json(work/'derivative_checks.json',checks)
                if not math.isfinite(error) or error>.01:raise ValueError('Alternative sampled Hessian column failed')
        record=dict(status='derivative_checked_solution_thermal',identity=dict(geometry_sha256=sha(out/'opt.xyz'),
            charge=0,mult=1,protocol=dict(SETTINGS,role='SVPD sensitivity; not the primary native protocol')),
            converged=True,stationary_minimum=True,grid_response=True,smd_spin_diagnostics=diag,
            gradient_max=float(abs(gradient).max()),maximum_asymmetry=asym,derivative_checks=checks,
            execution_memory=memory,software=installed_versions(),primary_native_certificate_claimed=False)
        atomic_write_json(native,record)
        _,partition=ensure_projected_thermal(native,out/'opt.xyz')
        if min(partition['frequencies_cm'])<-5:raise ValueError('Alternative internal minimum gate failed')
        primary=task['primary'];thermal=partition['g_thermal_total_eV'];total=energy['e_smd_eV']+thermal
        status.update(status='alternative_basis_control_checked',G_eV=total,E_eV=energy['e_smd_eV'],
            nuclear_by_cutoff_eV=partition['nuclear_by_cutoff_eV'],gradient_max=record['gradient_max'],
            minimum_internal_frequency_cm=min(partition['frequencies_cm']),
            delta_thermal_eV=thermal-primary['thermal_eV'],delta_geometry_SP_eV=energy['e_smd_eV']-primary['E_eV'],
            delta_first_reduction_V_holding_anion_and_Fc_fixed=total-primary['G_eV'],
            surface_record=str(native),partition_record=partition['postprocessing_path'],primary_native_admission=False)
        if any(sha(p)!=h for p,h in plan['source_sha256'].items()):raise ValueError('Control sources changed during calculation')
    except BaseException as error:
        status.update(status='failed_or_interrupted_preserving_records',error=repr(error));raise
    finally:
        status['finished']=time.time();atomic_write_json(out/'status.json',status)


def prepare(campaign):
    from redox.core.common import read_manifest
    from redox.core.protocol import load_record
    from redox.qm.native_thermal import load_native_thermal
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    rows={(r['id'],r['state']):r for r in read_manifest()};tasks=[]
    campaign.mkdir(parents=True,exist_ok=False)
    for key in PANEL:
        row=rows[key,'neu'];sd=DFT/key/'neu';raw=json.loads((sd/'result.json').read_text())
        energy=load_record(sd,0,1);thermal=load_native_thermal(sd,0,1)
        if energy is None or thermal is None:raise ValueError('Primary control state is not checked: '+key)
        target=campaign/(key+'.xyz');target.write_bytes((sd/'opt.xyz').read_bytes())
        tasks.append(dict(key=key,input=str(target),input_sha256=sha(target),smiles=row['smiles'],role='basis_method_control',
            primary=dict(E_eV=energy['e_smd_eV'],thermal_eV=thermal['g_thermal_total_eV'],G_eV=energy['e_smd_eV']+thermal['g_thermal_total_eV'],
                geometry_sha256=sha(sd/'opt.xyz'),state_dir=str(sd))))
    sources=[Path(__file__),*list((ROOT/'src/redox/qm').glob('*.py')),*list((ROOT/'src/redox/core').glob('*.py')),
        ROOT/'src/redox/workflows/geometry_checkpoints.py',ROOT/'config/projected_thermal_validation.json']
    atomic_write_json(campaign/'plan.json',dict(tasks=tasks,settings=SETTINGS,optimizer=OPTIMIZER,
        minimum_GPU_MiB=31744,raw_root=str(ROOT/'calcs/validation'/campaign.name),source_sha256={str(p):sha(p) for p in sources},
        comparisons='Neutral reoptimization plus SVPD Hessian; original anion and Fc fixed. Separate nuclear and geometry/SP effects. Not primary admission.'))
    print(json.dumps(dict(prepared=len(tasks),campaign=str(campaign))))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--prepare',action='store_true');parser.add_argument('--task');args=parser.parse_args()
    if args.prepare:prepare(args.campaign)
    elif args.task:run(args.campaign,args.task)
    else:parser.error('Use --prepare or --task')


if __name__=='__main__':main()
