# First admitted charged-polymer native component and matched followups

[Independent review](site2_resource_review.json) · [Seven actual starts](followup_start_verification.json)
· [Prepared followups](../jcim_completion_2026-11-02/native_followup_preparation.json)

The all-same, site-2 viologen component now passes independent native review:
114 atoms, global charge +2, singlet. Its complete Hessian, all twelve sampled
column/step comparisons, projected modes, symmetry, nuclear terms and component
sums reconstruct. The minimum internal frequency is +8.1893 cm⁻¹, maximum
sampled-column relative error 4.665e-6, and observed GPU peak 23,822 MiB. Normal
closure, actual SCF checkpoint/frame, gradients, grid-response, source bindings,
exclusive residency and the resource chain are checked separately from the
producer status. Its own G values are recorded by cutoff in the review.

This component uses r2SCAN-D4/def2-SVP electronic and Hessian terms under the
explicitly separate smooth-surface MeCN model. It is not a primary-model
ωB97M-V grafting/redox replacement, a converged ensemble or a balanced occupied
free-energy cycle. The review retains `primary_model_admission=false` and
`material_admission=false`. Absolute G values from different graphs are not
ranked as an interaction strength.

The passed pilot releases five other single-graft inputs and two matched empty
references whose dimensions and worker/kernel sources match its declared class.
Seven actual scientific processes were verified with PID/start identities,
disjoint two-CPU affinities, six thread caps=2 and one resident process per GPU:
five V100s on lambda2, one V100 on lambda4 and one H100 on lambda13. These are
running/unqualified followups; no completed G is inferred from the launch count.
The node with an earlier foreign-residency event was avoided for this pass.

An initial empty-reference placement was refused before scientific start by
aggregate RAM commitments. That controller was proven closed and its records
preserved. The H100 preflight then exposed three node-local library bindings
absent on that node. The [relocation](dependency_relocation.json) creates a fresh
namespace with byte-identical shared-storage replacements and changes no
geometry, state, method, source code or physical gates. A further reservation
path refusal revealed differing home-directory layouts; the common resource
layer now reads declared site defaults in `config/runtime_paths.json`. The
[successful H100 probe](lambda13_shared_probe.json) precedes its gated actual
start. These are dependency/prelaunch refusals, not failed physical calculations.

The original first-pilot namespace is
`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-01_viologen_all_same_occ_2_native_storage_recovery/`.
The six followup namespaces on V100 are under
`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-02_viologen_native_followup_prepared/`.
The shared-dependency H100 reference is at
`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-02_all_same_empty_shared_dependencies/`.
All previous attempts remain preserved.

All four double-graft components now pass stationary electronic review. The
all-same sites-2/4 case was recovered by a fresh final solve at its preserved
converged optimizer return; it did not repeat the optimization. The twelve
checked components give [four balanced electronic cycles](stationary_electronic_summary.json):
all-same adjacent −3.46 and one-spacer +18.68 kJ/mol; alternating adjacent
+1.20 and one-spacer −13.48 kJ/mol. These are stationary whole-chain E
nonadditivities under one model, not isolated motif binding or ΔG.

The 137-atom double grafts have 1,319 AOs and 249 occupied orbitals and exceed
the smaller pilot's dimension envelope. A separately declared H100 dimension
pilot retains the same compact solver, equations and physical thresholds, with
80-GiB-class device minimum, 88,000-MiB observed-peak ceiling and a 6-GiB retained
output envelope. Its submission/observation records are separate from completed
physical reviews. The old class remains capped at 31,744 MiB. The first preflight
refusal is preserved; new namespaces carry frozen review files to remain visible
inside the isolated worker checkout.

Minimum/Hessian and complete native-cycle gates remain mandatory after
stationarity. No spacing probability, graft loading, association free energy
or electrochemical transport is inferred.


[H100 actual-start verification](larger_native_actual_start.json) checks the
scientific PID/start, two disjoint CPUs, all six thread caps, sole GPU residency
and actual SCF progress. The [three prepared followups](larger_followup_preparation.json)
remain held for its physical/resource review. Their preparation is not a start.


The [eight-worker observation](eight_worker_observation.json) verifies all five
lambda2, one lambda4 and two lambda13 workers by live PID/start identity and
actual GPU residency, with disjoint CPU affinities. This is a dated observation,
not a completion certificate. The first component review's exact original
producer is retained under `producers/`, matching its reported SHA-256; the
current reviewer separately adds the larger resource-class checks.
