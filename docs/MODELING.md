# Modeling the system

How we turn optimized structures into stability + redox-potential descriptors. This is
the design we compute against; numbers/level-of-theory choices are flagged **[DECIDE]**.

## 1. States and the redox graph

Each molecule is a set of redox **states**, each with explicit `(charge, mult)`
(`config/redox_groups.py`, `library/manifest.csv`). A redox **event** connects two
adjacent states differing by one electron:

```
O + e-  ->  R        (reduction)      E° = -ΔG°/(1·F) - E_ref
```

Multi-electron groups are modeled as **sequential 1e events**, each with its own E°:
- viologen: `2+  --e-->  +•  --e-->  0`   → E°₁, E°₂
- anthraquinone: `0  --e-->  -1•  --e-->  -2`  → E°₁, E°₂

"Number of electrons transferred" (descriptor) = count of events inside an accessible
potential window.

## 2. Free energy of each state (solvated)

We optimize **in the SMD(MeCN) continuum**, so we get the solution-phase free energy
directly (no gas→solution thermodynamic cycle needed):

```
G_solv(state) = E_elec(SMD)  +  G_thermal
G_thermal     = ZPE + H_thermal(298K) - T·S          (from harmonic frequencies)
```

- `E_elec(SMD)`: DFT energy in the continuum at the SMD-optimized geometry.
- `G_thermal`: harmonic frequencies. **[DECIDE]** frequencies at opt level (rigorous) vs.
  a cheaper level (e.g. xtb) vs. ΔG_elec-only screening first (ignore thermal, add later).
  For charged organics the thermal term is often near-constant across a redox couple, so a
  first screen on E_elec alone is defensible; refine the shortlist with full G.

## 3. Redox potential

For event O + e⁻ → R:

```
ΔG°_solv = G_solv(R) - G_solv(O)          # electron free energy folded into referencing
E_abs    = -ΔG°_solv / F                  # absolute potential (V)
E_ref    = E_abs - E_ref(internal)        # report vs Fc/Fc+
```

**Referencing — compute the reference the same way.** Rather than trusting a literature
absolute SHE, we compute **ferrocene → ferrocenium** with the identical method/solvent and
reference every potential to it (Fc/Fc⁺). This cancels systematic method error and matches
common experimental practice in MeCN. Constants in `config/electrolyte.py`.

## 4. Reorganization energy λ (Marcus)

`λ = λ_inner + λ_outer`.

**Inner-sphere (4-point, from our optimizations):**
```
λ_i = [E_O(geom_R) - E_O(geom_O)] + [E_R(geom_O) - E_R(geom_R)]
```
Needs each state's energy at both geometries → 4 single points per event. **Cheap with
UMA** (charge/spin aware), which is where the MLIP earns its keep; cross-check on DFT.

**Outer-sphere (solvent):** continuum estimate (Marcus two-sphere, or nonequilibrium vs.
equilibrium SMD solvation). **[DECIDE]** two-sphere approximation vs. explicit
nonequilibrium-continuum from the DFT.

## 5. Structural descriptors

- **RMSD between redox states:** Kabsch-align optimized O vs R geometries (heavy atoms).
- **Structural change on ox/red:** key bond-length/angle deltas (e.g. C–N⁺, C=O, N–O),
  planarity/dihedral change of the redox core; report the largest movers.
- **Collapse check:** radius of gyration + intramolecular close-contact count vs. the
  ETKDG seed — flags artificial gas-phase folding (see build.py). Charged states are
  solvated-pre-optimized to prevent it.

## 6. Electronic descriptors

- HOMO / LUMO / gap; vertical & adiabatic IP/EA.
- Spin density localization (is the unpaired electron on the intended redox group?).
- Partial charges (e.g. Hirshfeld/CHELPG) and dipole; distributions across the library.

## 7. Functional-group stability

- **Connectivity preserved?** Re-perceive bonds from the optimized geometry; compare to
  input graph. A broken/rearranged group in some redox state = unstable in that state.
- **Energetic accessibility:** is the redox event in a sane potential window for MeCN?
- **Radical/ion stability:** spin contamination ⟨S²⟩, no spontaneous fragmentation,
  reasonable BDEs at the benzylic linker (the ex-Cl site is a known weak point).

## 8. Level of theory  **[DECIDE]**

Targets for comparability with **D3TaLES** (Risko/Odom, UKentucky) and with the UMA/OMol
reference (ωB97M-V/def2-TZVPD):

| Step | Candidate |
|------|-----------|
| Pre-opt (charged) | xtb GFN2 + ALPB(MeCN) |
| Pre-opt (neutral) | UMA (gas) |
| Opt + freq | ωB97X-D3(BJ) or B3LYP-D3 / def2-SVP(D) + SMD(MeCN) |
| Energy single point | ωB97M-V or ωB97X-D3 / def2-TZVP(D) + SMD(MeCN) |

Open decision: match D3TaLES's exact protocol for direct comparison, or use a
range-separated hybrid aligned with OMol. GPU acceleration via `gpu4pyscf` on the V100s.
