def dft_smd(xyz_path: Path, charge: int, mult: int,
            opt_xc: str = OPT_XC, opt_basis: str = OPT_BASIS,
            opt_basis_anion: str = OPT_BASIS_ANION, opt_disp: str | None = OPT_DISP,
            opt_nlc: str | None = OPT_NLC,
            sp_xc: str = SP_XC, sp_basis: str = SP_BASIS,
            sp_basis_anion: str = SP_BASIS_ANION, sp_disp: str | None = SP_DISP,
            sp_nlc: str | None = SP_NLC,
            solvent: str | None = None, do_gas: bool = True, do_smd: bool = True,
            do_opt: bool = True, do_freq: bool = True, opt_out: Path | None = None,
            max_opt_steps: int = 100, constraints: str | None = None,
            conv_params: dict | None = None, assert_convergence: bool = True,
            backend: str = "cpu") -> dict:
    """Composite DFT+SMD: optimize geometry in solvent at (opt_xc/opt_basis), then score
    the energy at (sp_xc/sp_basis). Anions use the diffuse basis variant. Returns energies
    (Ha) in SMD solvent and (optionally) gas phase. backend='gpu' uses gpu4pyscf.
    """
    from pyscf import gto
    from ase.io import read
    dft = _dft_module(backend)

    q = int(charge)
    nunpaired = int(mult) - 1  # pyscf mol.spin = 2S = (mult-1)
    obas = _basis_for(q, opt_basis, opt_basis_anion)
    sbas = _basis_for(q, sp_basis, sp_basis_anion)

    atoms = read(str(xyz_path))
    def _mol(basis):
        syms = atoms.get_chemical_symbols()
        atom_str = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                             for s, p in zip(syms, atoms.positions))
        return gto.M(atom=atom_str, basis=basis, ecp=_ecp_for(syms, basis),
                     charge=q, spin=nunpaired, verbose=0)

    solvent = solvent or _electrolyte().SOLVENT["name"]

    out = dict(xyz=str(xyz_path), charge=q, mult=int(mult), solvent=solvent,
               optimized=bool(do_opt), backend=backend,
               opt_xc=opt_xc, opt_basis=obas, opt_disp=opt_disp, opt_nlc=opt_nlc,
               sp_xc=sp_xc, sp_basis=sbas, sp_disp=sp_disp, sp_nlc=sp_nlc,
               # xc/basis mirror the ENERGY level (what redox referencing keys on)
               xc=sp_xc, basis=sbas)

    # Optimize the geometry IN SOLVENT (SMD gradients) at the opt level.
    mol = _mol(obas if do_opt else sbas)
    if do_opt:
        geomopt = _geomopt_fn(backend)
        mf_opt = _build_smd_mf(mol, opt_xc, solvent, backend, opt_disp, opt_nlc)
        # `constraints` (a geomeTRIC constraints file, e.g. frozen dihedrals) keeps a flexible
        # anion in the neutral's conformer -> a conformer-matched inner-sphere lambda.
        # `conv_params` (convergence_gmax/grms/energy...) + assert_convergence=False let a
        # constrained opt stop at a loose-but-energy-converged geometry instead of oscillating.
        _optkw = dict(maxsteps=max_opt_steps, assert_convergence=assert_convergence)
        if constraints:
            _optkw["constraints"] = constraints
        if conv_params:
            _optkw.update(conv_params)
        gradient = mf_opt.nuc_grad_method()
        gradient.grid_response = True
        scanner = gradient.as_scanner()
        converged, mol = geomopt(scanner, **_optkw)
        out.update(optimized=bool(converged and scanner.converged),
                   geometry_optimization_converged=bool(converged),
                   optimization_SCF_converged=bool(scanner.converged),
                   opt_gradient_grid_response=True)
        if assert_convergence and not out["optimized"]:
            raise RuntimeError("Geometry optimization did not converge; preserve the trajectory before resuming")
        opt_atoms = _mol_to_atoms(mol)
        if opt_out is not None:
            _write_xyz(opt_atoms, Path(opt_out),
                       comment=f"SMD({solvent})-opt {opt_xc}/{obas} q={q} m={mult}")
            out["opt_out"] = str(opt_out)
        # rebuild at the (larger) single-point basis on the optimized coordinates
        atoms = opt_atoms
        del gradient, scanner, mf_opt
        import gc
        gc.collect()
        _free_gpu_pool()

    # Energy: solvated single point at the sp level on the optimized geometry.
    # do_smd=False skips the (expensive) SMD SCF for callers that only need the gas energy
    # (e.g. inner-sphere reorganization cross-points) — a ~2x saving there.
    mol_sp = _mol(sbas)
    guesses = OPEN_SHELL_GUESSES if nunpaired else ("minao",)
    # provenance of the energy: what an adopter must match (redox.qm.sp._adoptable)
    from redox.core.protocol import software_versions
    from redox.core.protocol import geom_sha1
    scored_path = Path(opt_out) if do_opt and opt_out is not None else (Path(xyz_path) if not do_opt else None)
    out.update(density_fit=True, conv_tol=1e-9, max_cycle=200, scf_guesses=list(guesses),
               energy_protocol=dict(_SP, xc=sp_xc, basis=sbas, disp=sp_disp,
                                    nlc=sp_nlc, solvent=solvent),
               energy_geom_sha1=geom_sha1(scored_path) if scored_path is not None else None,
               software=software_versions())
    e_smd = None
    if do_smd:
        smd_seeds = (_metal_doublet_phase_seeds(mol_sp, sp_xc, backend, sp_disp, sp_nlc, solvent)
                     if nunpaired == 1 and 'Fe' in atoms.get_chemical_symbols() else {})
        e_smd, conv_smd, diag_smd, per = _lowest_scf(
            lambda: _build_smd_mf(mol_sp, sp_xc, solvent, backend, sp_disp, sp_nlc), guesses,
            post=lambda m: _orbital_diag(m, nunpaired), initial_densities=smd_seeds)
        out['smd_occupation_guesses'] = list(smd_seeds)
        out.update(e_smd_Ha=e_smd, e_smd_eV=e_smd * HARTREE_EV, converged_smd=conv_smd,
                   scf_guess_energies_smd=per)
        for k,v in (diag_smd or {}).items():
            out[f'smd_{k}']=v

    if do_gas:
        KS = dft.RKS if nunpaired == 0 else dft.UKS

        def _gas_mf():
            m = KS(mol_sp).density_fit(); m.xc = sp_xc   # RI-J: ~3-5x faster (see _build_smd_mf)
            if sp_disp:
                m.disp = sp_disp
            if sp_nlc:
                m.nlc = sp_nlc
            return m
        gas_seeds = (_metal_doublet_gas_seeds(mol_sp, sp_xc, backend, sp_disp, sp_nlc)
                     if nunpaired == 1 and 'Fe' in atoms.get_chemical_symbols() else {})
        e_gas, conv_gas, diag, per = _lowest_scf(_gas_mf, guesses,
                        post=lambda m: _orbital_diag(m, nunpaired), initial_densities=gas_seeds)
        out['gas_occupation_guesses'] = list(gas_seeds)
        out.update(e_gas_Ha=e_gas, e_gas_eV=e_gas * HARTREE_EV, converged_gas=conv_gas,
                   scf_guess_energies_gas=per)
        # QC: anion-bound check (gas HOMO > 0 => unbound) + spin contamination.
        for k, v in (diag or {}).items():
            out[f"gas_{k}"] = v
        if out.get("gas_homo_eV") is not None and q < 0:
            out["anion_unbound"] = bool(out["gas_homo_eV"] > 0.0)
        if e_smd is not None:
            out["dG_solv_eV"] = (e_smd - e_gas) * HARTREE_EV

    # Thermal free-energy correction (GFN2-xTB RRHO on the DFT+SMD-optimized geometry) so
    # redox uses G = E_smd + g_thermal, not the bare electronic energy. See _thermal_correction.
    if do_freq:
        try:
            th = _thermal_correction(atoms, q, int(mult))
            out.update(**th, freq_level="gfn2-xtb (RRHO, 298.15K)")
        except Exception as exc:  # a failed Hessian must not lose the (expensive) energies
            print(f"[warn] freq/thermal failed: {type(exc).__name__}: {str(exc)[:120]}",
                  flush=True)
            out["g_thermal_eV"] = None
    return out
