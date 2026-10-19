# Parallel validation and candidate expansion — 2026-10-10

Two efforts are running concurrently. Candidate expansion has produced a reviewed proposal
pool; the bounded quantum repair campaign is still computing. The first preserved table
snapshot passed 8,306 independent checks. It is not the final post-repair calculation set.
Current code passes **58 tests**.

## 🔄 Continuation on 2026-10-11

The current audited checkpoint is
`results/review_2026-10-11_spin_and_thermal_checkpoint/`: **10,223 checks, zero discrepancies**.
`results/LATEST_REVIEW.json` identifies the current snapshot; older root-level exports are
preserved. Imaginary-mode RRHO contributions now appear only in explicitly named potential
and disproportionation **estimate** columns. Their ordinary property columns are blank and
their rows are `THERMAL_UNVALIDATED`; incomplete paths never enter the primary Pareto front.
The original estimates remain available for diagnosis.

Completed physical repairs and new evidence:

- The incorrect Fc+ gas own and cross energies were lowered by 2.536 and 2.354 eV.
  An occupation-seeded production policy reproduces the independent repaired own energy
  within 2.2e-11 eV, preventing recurrence. SMD Fc referencing is unchanged.
- `d3_80izcq/red2` was optimized in the wrong spin state: the singlet is **1.757 eV below**
  the UMA-selected triplet at active SMD energies. The singlet is now chosen, all old files
  are preserved under `physical_thermal_spin_2026-10-10/superseded_spin_states/`, and a
  species/protocol-checked override prevents future reruns from restoring the triplet.
  Its old one-electron capacity classification is not trusted; second-wave G is incomplete
  pending thermal validation. The new state's thermal correction was not copied from the
  triplet. Old-spin cross caches do not match its new geometry/spin hashes.
- Electronic singlet-triplet checks confirm `aq_benzylamino`, `aq_benzyloxy`, and PMDI
  singlets with margins 0.795, 0.864 and 0.720 eV. Remaining sub-0.5 eV UMA cases are in
  a bounded one-GPU queue. Small unverified UMA gaps and known lower alternate states
  cannot qualify a row for the complete primary front.
- Three nonequilibrium PCM controls completed with independent contraction agreement.
  Symmetric one-body electronic solution relaxations are 1.152 eV for parent AQ's first
  reduction, 1.121 eV for its dianion step, and 1.146 eV for the grafted AQ first reduction.
  These are distinct from gas inner-only lambda and remain diagnostic pending broader
  validation; they are never added to the old outer term.
- Grafted AQ's neutral-like local red1 conformer is 0.0232 eV above its original red1
  geometry; a linker-fixed version is 0.0240 eV above. These are explicit local/fixed-linker
  sensitivity tests, not proof that a polymer freezes the linker or that flags should vanish.
- Gas GFN2 produces two deep AQ-dianion imaginary modes localized on the tether. ALPB
  at the same coordinates removes those deep modes. Mode-following in solution found
  genuine GFN2-ALPB minima for the tested AQ neutral and dianion. Transfer accuracy of
  low-level thermal terms to the high-level energy still requires the DFT reference.
- The installed analytic DFT Hessian failed two-size gradient tests by up to 27% and was
  rejected. A full numerical Hessian was completed, and a small nonconservative component
  was traced to omitted integration-grid response in gradients. Including that response
  reduces the worst-pair asymmetry from 1.78e-3 to 8.46e-9 Ha/Bohr²; level-4 grids give
  consistent checked columns. The full conservative-gradient Hessian is being recomputed
  with preserved data, independent cold-start/step checks and two reserved GPU stages.
- A source-level OROP check corrects Fe's LANL2DZ basis/ECP and distinguishes electronic
  from legacy RRHO reference values. Fc is 4.361 V electronically under the partial
  reproduced protocol versus published 4.662 V. Basis/ECP and thermal convention account
  for approximately 0.18 V of the earlier 0.48 V discrepancy. The residual ~0.30 V is
  unresolved; no fitted correction or external Fc constant was substituted.

The two new registered pilot motifs, ether viologen and N-benzyl phthalimide control,
have **five converged UMA charge states** with composition, connectivity and spin-parity
checks passed. Their results are gas warm starts and not solution property predictions.
The full six-motif, 17-state isolated preparation remains available. The newly documented
NDI terminal-alcohol precursor offers stronger synthetic evidence than speculative mono-N-H
NDI grafting; its proposed ether transfer to resin chloride remains to be established.

Live status and evidence:

- `results/validation/physical_resume_2026-10-11_manifest.json`
- `results/validation/physical_thermal_spin_2026-10-10/all_small_spin_checks/worker_status.json`
- `results/validation/physical_thermal_spin_2026-10-10/numerical_solution_hessian_gridresponse/`
- `results/validation/physical_thermal_spin_2026-10-10/numerical_solution_hessian_gridresponse_verified/`
- `results/discovery/uma_pilot_2026-10-11/status.json` and `geometry_checks.json`
- `results/discovery/next_batch_preparation_2026-10-10/BATCH_READINESS.md`

Completion requires the conservative numerical Hessian checks, remaining spin checks,
new-state thermal/cross calculations, and a fresh independent export audit. Arithmetic
agreement alone does not establish experiment-level physical accuracy.

## 🔬 Metric validation

Ten missing active-protocol energy records were recovered from existing `result.json`
calculations that prove the requested protocol and SCF provenance. Their source geometry
and result-file SHA256 hashes were unchanged. One previously declared state,
`d3_05hfap/red1`, lacks a completed DFT result because earlier 15 GB GPU runs exhausted
memory after optimization. It is being repaired on a reserved 32 GB V100.

| Worker | Host / GPU | CPU affinity | Authorized work | State at this update |
|---|---|---|---|---|
| Geometry | lambda1 / 0 | 0–1 | Three-seed optimization, thermal term, active SP and four cross points for `d3_05hfap` | First seed optimized; production SP running |
| Cross points | lambda1 / 1 | 2–3 | 24 missing gas cross points for six existing discovery molecules | First record verified; computation running |
| Finalizer | lambda1 / none | 4 | Wait, verify eight discovery molecules, export and independently audit new tables | Waiting for workers; six-hour deadline |

The finalizer adds no quantum workload and does not retry failures. It records failures
or timeout instead of claiming completion. It will write a fresh snapshot to
`results/validation/parallel_metrics_2026-10-10/after_repair_tables` only after coverage and
cross-cache checks pass; completion also requires the independent table/workbook audit.

The physical audit established correct atomic composition and spin parity for all 196
existing geometries. All 392 gas/SMD phase energies have consistent Hartree/eV units;
282 phases explicitly equal the lowest recorded converged SCF guess. The remaining 110
have no retained per-guess evidence and cannot be verified on that point from the JSON.
There are 38 gas-unbound states and 16 imaginary-mode thermal states; these remain explicit
scientific limitations. The small UMA spin gap for `aq_benzylamino/red2` still needs a DFT
singlet/triplet check.

The methyl-viologen/PF6 helper `mv_ip1/s0` has changed chemical connectivity: the optimized
structure contains new C–F/C–P bonds and lost P–F connectivity. Comparison with its original
library seed confirms that the intended ion pair was not retained. This helper was already
excluded from candidate scorecards. Readiness now labels it **INVALID**, rather than allowing
converged energies to imply a chemically correct species.

The post-adoption readiness scan with structure checks reports 179 ready, 16 REVIEW,
one INCOMPLETE and one INVALID out of 197 declared states. This is a checkpoint while the
missing geometry is being repaired, not a prediction of the final counts.

Additional software fixes require finite cross energies and full embedded protocol
provenance, inspect every anionic leg of the four-point calculation, strip radical labels
when matching geometry skeletons, and map constrained torsions into the actual XYZ atom
order. Equilibrium-SMD conformer distortion energies remain diagnostic and are not chosen
as a validated inner-sphere reorganization value. Thermal sensitivity caches now depend on
geometry/charge/spin/cutoff, reject failed xTB runs, and use Fc at the same cutoff; no new
thermal sensitivity calculation was submitted during this campaign.

Evidence and live status:

- [`launch_manifest.json`](../results/validation/parallel_metrics_2026-10-10/launch_manifest.json)
- [`finalization_status.json`](../results/validation/parallel_metrics_2026-10-10/finalization_status.json)
- [`REPAIR_AND_PHYSICAL_REVIEW.md`](../results/validation/parallel_metrics_2026-10-10/REPAIR_AND_PHYSICAL_REVIEW.md)
- [`readiness_with_structure.csv`](../results/validation/parallel_metrics_2026-10-10/readiness_with_structure.csv)
- Logs under `logs/validation_parallel_2026-10-10/`.

## 🧪 Candidate expansion

The refreshed D3TaLES pool has 104 new-to-library quinone molecules. Further chemistry and
literature review produced **39 proposals, representing 38 connectivities**: one PMDI
connectivity has explicit R/S alternatives. This includes 20 additions to the earlier
19-molecule preliminary set.

| Family | Proposals | Role |
|---|---:|---|
| Quinones | 19 | Anolyte candidates with attachment and competing-reaction flags |
| Diimides | 15 | PMDI/NDI side-chain and attachment controls |
| Ether viologens | 2 | Two-electron anolyte designs; PF6 performance unverified |
| Monoimides | 2 | One-electron controls |
| Dimethoxyphenothiazine | 1 | Specific source-inspired two-oxidation catholyte proposal |

All 39 pass state/electron parity checks, fully specified stereochemistry, precursor-to-graft
round trip, carbonyl preservation, real ETKDG embedding and the builder dry-run. An independent
atom-sum mass/capacity audit passes 156 checks. These tests establish consistent proposed
structures, not experimental synthetic feasibility or computed production properties.
Parent literature evidence is separated from the uncomputed grafted derivatives.

The first eight proposals for chemistry/compute staging are:

1. `design_pmdi_methoxyethyl`
2. `design_ndi_methyl`
3. `design_ndi_methoxyethyl`
4. `d3_80prff`
5. `d3_80eeoz`
6. `exp_aq_1oxy`
7. `exp_aq_2oxy_6methoxy`
8. `exp_viologen_methoxyethyl`

The ten earlier D3TaLES front entries require route review: a substructure reaction does
not prove selective Merrifield functionalization. Conjugated aminoquinone nitrogens may be
weak nucleophiles, and aza-quinones can compete for ring-N alkylation. Ether-capped imides
are a controlled polarity/mass/conformer series, not an assumed solubility improvement.
The substituted phenothiazine branch is supported by specific two-electron parent work; it
does not change the one-electron assumption for unsubstituted phenothiazine.[^1]

No new proposal has yet been registered or submitted for production QM. Existing-cohort
metric repair and proposal expansion are independent; all future production candidates
must satisfy the same structure, energy, integrity and ranking checks.

Artifacts:

- [`CANDIDATE_REVIEW.md`](../results/discovery/parallel_expansion_2026-10-10/CANDIDATE_REVIEW.md)
- [`curated_candidates.csv`](../results/discovery/parallel_expansion_2026-10-10/curated_candidates.csv)
- [`proposals.py`](../results/discovery/parallel_expansion_2026-10-10/proposals.py)
- [`sources.json`](../results/discovery/parallel_expansion_2026-10-10/sources.json)
- [`validation_summary.json`](../results/discovery/parallel_expansion_2026-10-10/validation_summary.json)

[^1]: [A stable two-electron-donating phenothiazine for application in nonaqueous redox flow batteries](https://pubs.rsc.org/en/content/articlelanding/2017/ta/c7ta05883g); [Tailoring Two-Electron-Donating Phenothiazines To Enable High-Concentration Redox Electrolytes](https://pubs.acs.org/doi/10.1021/acs.chemmater.8b04770). Exact proposed grafted structures have no verified experimental performance.
