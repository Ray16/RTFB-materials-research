def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mult", type=int, choices=(1, 2), required=True)
    parser.add_argument("--workspace-GiB", type=float, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    verify_gate_owner()
    caps = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(key) != "2" for key in caps):
        raise ValueError("Control must inherit two disjoint CPU cores and matching native limits")
    import pyscf.lib
    from pyscf import gto
    from redox.qm import dft
    from redox.qm.cds_finite import finer_cds_hessian
    pyscf.lib.num_threads(2)
    source = ROOT / "scripts/pipeline/df_occupied_workspace.py"
    spec = importlib.util.spec_from_file_location("qualified_SMD_buffer_adapter", source)
    memory = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(memory)
    args.out.mkdir(parents=True, exist_ok=False)
    xyz = args.out / "input.xyz"
    xyz.write_text("4\nSinglet implementation control\nC 0 0 0\nO 1.21 0 0\nH -.61 .94 0\nH -.61 -.94 0\n" if args.mult == 1 else "4\nDoublet implementation control\nC 0 0 0\nH 1.08 0 0\nH -.54 .935307436 0\nH -.54 -.935307436 0\n")
    report = dict(status="comparing_original_and_compact_execution", mult=args.mult, host=os.uname().nodename,
        pid=os.getpid(), process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), GPU=os.environ["CUDA_VISIBLE_DEVICES"], started=time.time(), arms={})
    atomic_write_json(args.out / "status.json", report)
    try:
        for name in ("original", "compact"):
            arm = args.out / name
            arm.mkdir()
            context = memory.compact_smd_execution(arm / "execution", budget_GiB=args.workspace_GiB, occupied_cap=4) if name == "compact" else nullcontext()
            with context:
                energy = dft.dft_smd(xyz, 0, args.mult, do_opt=True, do_freq=False, backend="gpu",
                    opt_out=arm / "opt.xyz", max_opt_steps=100)
                if not (energy["geometry_optimization_converged"] and energy["optimization_SCF_converged"] and energy["converged_smd"] and energy["converged_gas"]):
                    raise ValueError("Control optimization/SP did not converge")
                atoms = read(arm / "opt.xyz")
                mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(), atoms.positions)), basis="def2-svp", charge=0,
                    spin=args.mult - 1, unit="Angstrom", verbose=0, max_memory=6000)
                mf = dft._build_smd_mf(mol, "r2scan", "acetonitrile", "gpu", disp="d4")
                mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-12, 1e-8, 200
                mf.grids.level = 3
                mf.kernel()
                if not mf.converged:
                    raise ValueError("Geometry-level control SCF did not converge")
                derivative = mf.nuc_grad_method()
                derivative.grid_response = True
                gradient = derivative.kernel()
                gradient = np.asarray(gradient.get() if hasattr(gradient, "get") else gradient)
                mf.cphf_grids = mf.grids.copy()
                hessian = mf.Hessian()
                hessian.grid_response = True
                hessian.max_cycle = 100
                with finer_cds_hessian(1e-5):
                    matrix = hessian.kernel()
                matrix = np.asarray(matrix.get() if hasattr(matrix, "get") else matrix)
                np.save(arm / "gradient.npy", gradient)
                np.save(arm / "hessian.npy", matrix)
                atomic_write_json(arm / "energy.json", energy)
                report["arms"][name] = dict(E_smd_eV=energy["e_smd_eV"], E_gas_eV=energy["e_gas_eV"],
                    geometry_level_E_Ha=float(mf.e_tot), coordinates_A=atoms.positions.tolist())
                del hessian, derivative, mf, mol
                dft._free_gpu_pool()
            atomic_write_json(args.out / "status.json", report)
        original, compact = report["arms"]["original"], report["arms"]["compact"]
        differences = dict(
            maximum_energy_difference_eV=max(abs(original[key] - compact[key]) for key in ("E_smd_eV", "E_gas_eV")),
            geometry_level_energy_difference_Ha=abs(original["geometry_level_E_Ha"] - compact["geometry_level_E_Ha"]),
            maximum_coordinate_difference_A=float(np.max(abs(np.array(original["coordinates_A"]) - compact["coordinates_A"]))),
            maximum_gradient_difference=float(np.max(abs(np.load(args.out / "original/gradient.npy") - np.load(args.out / "compact/gradient.npy")))),
            maximum_Hessian_difference=float(np.max(abs(np.load(args.out / "original/hessian.npy") - np.load(args.out / "compact/hessian.npy")))),
        )
        if any(not np.isfinite(value) or value > 1e-6 for value in differences.values()):
            raise ValueError("Full energy/geometry/gradient/Hessian equivalence failed")
        import hashlib
        sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
        report.update(status="full_original_compact_equivalence_passed", differences=differences,
            evidence_sha256={str(p): sha(p) for p in args.out.rglob("*") if p.is_file() and p.name != "status.json"},
            source_sha256={str(p): sha(p) for p in (Path(__file__).resolve(), source, Path(dft.__file__))},
            scientific_model_changed=False)
    except BaseException as error:
        report.update(status="control_failed_outputs_preserved", error=repr(error))
        raise
    finally:
        report["finished"] = time.time()
        atomic_write_json(args.out / "status.json", report)
