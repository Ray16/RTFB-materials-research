# Current results

Use [INDEX.json](INDEX.json) for machine-readable current pointers. Historical
scientific outputs remain immutable at their original locations, including
`/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/`. Startup records alone do not
establish that a calculation is currently running.

- [Paper and SI](../manuscript/README.md): checked 21/30-page working documents,
  frozen evidence, independent arithmetic and a self-contained upload bundle.
  The current revision integrates both empty-backbone sampling comparisons and
  portable numerical replay of its 75 native states, plus six fresh checked
  TEMPO solution-response calculations and the checked exact viologen ladder.
- [Exact-molecule CPU example](../examples/exact_viologen_replay/README.md):
  ten versioned cached states reconstruct both viologen reductions and local
  quaternization with NumPy only. Nine Hessians and six saved-gradient panels
  pass extracted-archive replay with original scientific paths blocked.
- [Focused slides](../slides/README.md): fifteen main slides and four supporting
  slides, including the UMA-to-DFT pipeline, grafting ΔG and exact TEMPO response.
- [Exact TEMPO follow-up](../validation/tempo_solution_pair_2026-10-30/README.md):
  checked MeCN relaxation costs of 1.477/1.639 eV. The 10.4% directional asymmetry
  leaves a symmetric parameter unqualified; the all-axis material decision is pending.
- [Candidate pool](../validation/candidate_pool_2026-10-29/README.md): 904 exact
  local-graft graphs, with source/recipe identity and prospective route status.
  Two exact graphs now have reviewed route passes (TEMPO and baseline viologen);
  902 remain on hold. The [new exact-route review](../validation/viologen_route_2026-10-30/README.md)
  preserves the prior registry and all graph/property records.
  Enumeration and recorded property presence are separate from quantum or
  experimental validation; no complete material shortlist is claimed.
- [Exact enumerated viologen](../validation/viologen_native_recovery_2026-10-30/README.md):
  three native states pass independent review after minimum repair. MeCN reductions
  are −0.608/−1.208 V versus Fc; modeled monomer disproportionation is +57.91 kJ/mol.
  The [new quaternization cycle](../validation/viologen_grafting_2026-10-31/README.md) gives −11.49 kJ/mol in MeCN at 298.15 K from three checked precursor seeds. This is a finite local-site separated-ion model; polymer occupancy/coupling and material admission remain incomplete.
- [Charged viologen polymer preparation](../validation/viologen_polymer_preparation_2026-10-31/README.md):
  twelve exact graph/configuration models pass initial XYZ connectivity,
  stereochemistry, backbone and salt-mass audits. The
  [construction figure](figures/polymer_models/viologen_occupancy_construction.pdf)
  shows adjacent and one-spacer grafts. Charged-polymer thermochemistry and
  converged conformer ensembles remain incomplete; the exact-candidate material decision remains pending.
- [Charged viologen geometry campaign](../validation/viologen_polymer_geometry_2026-10-31/README.md):
  ten actual SCF starts verified across ten V100 GPUs, covering three single-site
  and adjacent/spacer graft patterns in each of two backbone configurations.
  Each optimizer retains evaluated geometries and a rolling SCF checkpoint.
  The independent final-geometry reviewer and 62 focused tests are checked;
  optimized structures, Hessian minima and interaction/free-energy conclusions
  remain pending. The frozen start observation is separate from live progress.
- [Grafting and candidate evidence](../validation/candidate_grafting_audit_2026-10-28/README.md):
  source-bound PMDI arithmetic, defined routes and explicit limits. The
  [DBUH–chloride check](../validation/dbu_chloride_anchor_2026-10-29/README.md)
  remains preliminary; its experimental conventions are not fully resolved.
- [Empty-backbone sampling](../validation/polymer_empty_batches_2026-10-30/README.md):
  forty native inputs in two configurations yield 20/19 distinct structures.
  Nominal 10-to-20 ensemble changes are −1.03/−1.23 kJ/mol; both empty components
  meet the declared finite-batch criterion. The linked native census contains
  75 reconstructed states. Occupied ensemble cycles remain incomplete.
- [Screening implementation assessment](../validation/screening_chemical_impact_2026-10-31/README.md):
  the installed GPU4PySCF screening defect reproduces on V100/H100; the isolated
  upstream correction passes both integration controls and V100/H100 closed/open-shell
  MeCN derivative controls. All 71 benchmark electronic states pass independent
  repeats, with all 39 electronic relative-potential shifts below 9e-11 V.
  Fresh and correctly scoped full polymer-Hessian restart controls are active.
  Complete polymer/reaction-cycle and production qualification remain open; new uncorrected admissions
  are held and stored chemical results remain provisional for release.
- [Working occupied-state census](../validation/screening_implementation_2026-10-31/native_census.json):
  122 input states now reconstruct under their stored protocol. Six earlier states verify
  14,643 source bindings, including one recovered interrupted case. The
  [coverage figure](figures/polymer_models/2026-10-31_response_recovery/native_occupancy_coverage.pdf)
  shows all twelve graph/configuration cells. Occupied ensemble cycles remain
  incomplete; the article retains its separate 75-state numerical snapshot.
  At observed class maxima, unfinished inputs have a full retained footprint of
  about 100 GiB, including some already written partial files. Shared storage has
  remained below the legacy new-start floor at the latest campaign observation;
  the measured sizes are not peak bounds. The bounded charged-geometry class has
  separate prospective output budgets and does not relax that legacy gate.
- [V100 response-execution control](../validation/screening_chemical_impact_2026-10-31/response_scope.json):
  a repeatable CPU allocation test supports disabling NumPy huge-page advice on
  the affected host. Three traced response probes failed and remain unqualified;
  the earlier seven-block/315-coordinate calculation used solute H1 and frozen
  solvent response. Its file/frame/resource checks are instrumentation evidence;
  they do not qualify a full equilibrium-solvent response or native G.
- [Response-checkpoint continuation](../validation/response_continuation_2026-10-31/README.md):
  full-solvent source/frame and saved-prefix checks pass for both live controls;
  98 focused tests check continuation, residuals, precision and ownership.
  The new continuation path requires a closed source and a fresh unchanged-equation
  residual check. Molecular continuation equivalence and full-Hessian qualification
  remain pending; no source worker was interrupted or new GPU job submitted.
  The [prospective full-Hessian comparison](../docs/HESSIAN_REPRODUCTION_CHECK.md)
  fixes numerical reproduction limits before the completed matrices are available.
- [Portable native numerical cache](../validation/native_cache_2026-10-30/README.md):
  a 63 MB raw-input archive for the same 75 states passes extracted-archive
  replay in a fresh hash-locked NumPy environment with original paths blocked.
  Complete historical runtime closure and public deposition remain open.
- [Remaining occupied-energy workflow](../validation/polymer_energy_classes_2026-10-30/README.md):
  all 160 declared inputs are represented by 156 prepared jobs and four identical
  checked references. Consult the linked live dispatcher and independent
  collector for current evaluation counts. Energy completion does not supply
  missing native thermochemistry.
- [Checked-energy to native thermochemistry](../validation/native_conveyor_2026-10-30/README.md):
  automatic preparation, guarded dispatch and independent reconstruction of the
  remaining occupied states. Ready thermochemistry has admission priority; the
  [40-GiB pilot](../validation/polymer_effective_checked_2026-10-30/README.md)
  now passes both physical and resource checks for one occupied conformation.
  Wider matched-class use still requires prospective checks.
- [Matched native resource profiles](../validation/native_profile_dispatch_2026-10-30/README.md):
  seven pending profiles prepared and two new V100 starts checked, including actual
  SCF/gradient checkpoints. The scheduling update retains the original numerical
  gates and adds shared-volume headroom before new starts.
- [Resource-refusal recovery](../validation/resource_refusals_2026-10-30/README.md):
  native and electronic prelaunch refusals now recover after source/ownership
  checks and a full policy-node producer audit, without charging scientific attempts.
- [Route and recent polymer literature](../validation/route_evidence_followup_2026-10-30/README.md):
  two close 2025/2026 studies are integrated into the paper; exact-route holds
  and numerical-model limits remain explicit.
- [Retained-state physical replay](../validation/publication_native_replay_2026-10-28/README.md)
  and [subsequent checked continuations](../validation/polymer_continuations_checked_2026-10-30/README.md)
  preserve the source bindings underlying the current collection.
- [Exact-SI energy conversion](../validation/energy_conversion_2026-10-30/README.md)
  retains historical arithmetic explicitly. The
  [saved-checkpoint counterpoise diagnostic](../validation/counterpoise_checkpoint_diagnostic_2026-10-30/README.md)
  remains separate from a completed five-term interaction result.

[Publication acceptance criteria](../docs/PUBLICATION_READINESS.md) and
[property-validation priorities](../docs/PROPERTY_VALIDATION_AND_PUBLICATION.md)
state the remaining scientific and public-release requirements. The current
Nitro force-repeatability failure remains explicit; gas-phase structural warm
starts do not replace final DFT/acetonitrile energies or Hessians.
