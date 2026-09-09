# Modeling the system

How we turn optimized structures into stability + redox-potential descriptors. This is
the design we compute against; numbers/level-of-theory choices are flagged **[DECIDE]**.

## 0. Candidate set: standalone (before) vs grafted (after functionalization)

The screening targets **polymer-bound** redox-active materials, so every candidate is modeled
as a **grafted model-compound**: the redox core carries its real functional substituent on one
side and a **4-methylbenzyl group as a proxy for the tether to the monomer/backbone** (the
"Merrifield handle") on the other. All six candidates are modeled this way, uniformly:

| Candidate id | Class | Grafted form (modeled) | Standalone id (before functionalization) |
|---|---|---|---|
| `viologen` | viologen | N-(4-methylbenzyl)-N′-**methyl** bipyridinium | `methyl_viologen` (N,N′-dimethyl) |
| `ethylviologen` | viologen | N-(4-methylbenzyl)-N′-**ethyl** bipyridinium | `ethylviologen_sa` (N,N′-diethyl) |
| `pmdi` | imide | N-(4-methylbenzyl)-N′-methyl pyromellitic diimide | `pmdi_sa` (N,N′-dimethyl) |
| `ndi_ammonium` | imide | N-(4-methylbenzyl)-N′-(ammoniopropyl) NDI | `ndi_ammonium_sa` |
| `mophquinone` | quinone | 2-(4-benzyloxyphenyl)-1,4-benzoquinone | `mophquinone_sa` (2-(4-methoxyphenyl)-, = xlsx/D3TaLES) |
| `dmophquinone` | quinone | 2-(5-benzyloxy-2-methoxyphenyl)-5-methoxy-quinone | `dmophquinone_sa` (all-methoxy, = xlsx/D3TaLES) |

- **"Before functionalization" = standalone (`*_sa`, and `methyl_viologen`)** — the bare redox
  core with a minimal cap (methyl, or the native methoxy for the quinones). These are NOT
  screening candidates; they are the reference row of the standalone-vs-grafted 2×2 and double
  as validation refs (the two quinone standalones are *exactly* the Candidates.xlsx / D3TaLES
  molecules, enabling a direct reorg cross-check). Tagged `family="validation"` in the manifest.
- **"After functionalization" = grafted** — the six ids above; these carry the SA/λ numbers
  reported in the candidate figure and scorecard.

**Naming convention in figures:** short class names ("methyl viologen", "ethyl viologen",
"PMDI", "ammonium-NDI", "MeO-phenyl quinone", "(MeO)₂-phenyl quinone") label the *grafted*
candidates; the grafted-model + 4-methylbenzyl-proxy convention is stated once in the caption.
`Candidates.xlsx` gives these as neutral, mono-substituted schematics — the SMILES there encode
the *substituent* (methyl / ethyl / methoxy), and the real modeled species is the charged,
fully-substituted structure above.

Why this is reasonable: the redox core dominates the inner-sphere λ, so a small cap faithfully
reproduces the grafted unit's reorganization energy at low cost; the benzyl proxy best matches
an aromatic/styrenic backbone. Caveats: grafted models are **asymmetric** (the two 1e steps
become slightly inequivalent vs. the symmetric standalone), and the benzyl adds a rotatable
bond → a real (small) conformational contribution to λ (see the torsion-scan check; grafting
barely moves methyl viologen's λ, ±5 meV, but shifts ethyl viologen's λ(+1/0) by ~+54 meV).

## 1. States and the redox graph

Each molecule is a set of redox **states**, each with explicit `(charge, mult)`
(`config/redox_groups.py`, `library/manifest.csv`). A redox **event** connects two
adjacent states differing by one electron:

```
O + e-  ->  R        (reduction)      E° = -ΔG°/(1·F) - E_ref
```

**Spin state is DETERMINED, not assumed.** For each charge we scan candidate multiplicities
with UMA (even-electron → singlet+triplet; odd → doublet; plus the config hint) and take
the lowest-energy one. The singlet–triplet gap is recorded (`spin_gap_eV`). This matters for
even-electron reduced species (viologen 2e-reduced, anthraquinone dianion), which could be
closed-/open-shell singlet or triplet — the config `mult_hint` is only a starting guess.

**Conformers are sampled per redox state.** build.py emits a conformer ENSEMBLE (top-K,
0.1 Å pruning); UMA relaxes every (conformer × multiplicity) for each state and keeps the
global minimum (`conf_idx`). FF ranking is only a coarse pre-filter — UMA (charge/spin-aware)
does the real ranking, since the preferred conformation changes between oxidation states.

**Protonation state — assumption (state explicitly).** We model **aprotic acetonitrile,
pure outer-sphere electron transfer** with **fixed protonation** (no protons added/removed
between redox states). This is appropriate for dry MeCN. CAVEAT: **proton-coupled electron
transfer (PCET)** can dominate for **quinone radical anions/dianions** (→ semiquinone/
hydroquinone) and **TEMPO reduction** (→ hydroxylamine) if any proton source (trace water,
acidic electrolyte) is present; then the measured potential differs from our pure-ET value.
Treat those groups' reduced states with this caveat; a separate PCET calculation is needed
if the electrolyte is not rigorously dry.

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

## 8. Higher-order effects & limitations

**π–π interactions / radical dimerization.** Several groups have aromatic redox states
that stack: viologen⁺• (classic π-dimer), phenothiazine⁺•, quinone radicals. Two regimes:
- *Intramolecular* (benzyl handle folding onto the redox aromatic) — captured by the
  single-molecule model; watched by the collapse guard (Rg + contact count, §5).
- *Intermolecular* (stacking/dimerization between neighboring groups in the packed resin)
  — **NOT captured** by the single-molecule + implicit-solvent model. This can shift the
  effective E° (merge/split multi-electron waves, drive comproportionation), especially
  for the dimerization-prone radicals above.

Implications:
- **Requires a dispersion-corrected functional** (ωB97M-V/VV10 or ωB97X-D); plain B3LYP
  misses π–π. Reinforces the §8 level-of-theory choice.
- Reported single-molecule E° is the *isolated-monomer* value; note aggregation as a
  caveat in P5. Optional follow-up: explicit π-dimer calcs (UMA-screened) for viologen⁺•
  / phenothiazine⁺• → dimerization energy + ΔE° shift, checked vs the known viologen case.

## 9. Level of theory  **[DECIDE]**

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
