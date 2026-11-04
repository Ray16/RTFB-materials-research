# Phenothiazine: oxidized stationarity failure and strict refinement

All three strict geometry refinements now pass independent review, including the
[neutral plus branch](../phenothiazine_oxidized_geometry_2026-11-02/neutral_plus_geometry_review.json).
The [fresh electronic campaign](../refined_electronic_2026-11-02/README.md) retains
both neutral branches and the oxidized state; one actual V100 start is verified.
Own electronic/native qualification remains pending. Earlier observations below
retain their original waiting/running counts.

The [oxidized repair](oxidized_geometry_review.json) and
[neutral minus branch](neutral_minus_geometry_review.json) now pass independent
strict geometry review of all 43 and 53 saved SCF frames, respectively. Their
fresh final gradients are 1.05e-6 and 3.50e-7 Ha/Bohr, and peak GPU memory is
9,118 and 8,632 MiB. Exact frames, graph/hydrogen topology, occupied spin,
method and closed resource/ownership records pass. The
[second neutral branch](second_neutral_actual_start.json) has a verified actual
start on lambda1. These are geometry qualifications: new high-level electronic
energies, full own Hessians/minimum checks and sulfur-column validation remain
required before G or a potential. Both neutral branches remain in the campaign.

The [independent diagnosis](stationarity_diagnosis.json) reconstructs the actual
oxidized +1 doublet SCF, wavefunction, exact geometry and gradient. The maximum
gradient is 0.00045735997 Ha/Bohr, above the unchanged 0.00045 limit. This attempt
stopped before a Hessian, sulfur-column validation or nuclear partition. Its
closed scientific, wrapper and budget owners, exclusive GPU residency and host
resource chain distinguish this physical failure from the earlier interrupted
CUDA startup. Neither attempt supplies G or an oxidation potential.

The [current record](current.json) links a fresh, strict geometry-only refinement
of that exact 39-atom/397-AO +1 doublet in primary SMD MeCN. Corrected
r2SCAN-D4/def2-SVP, the original tight optimizer thresholds and validated 32-GB
device class are retained. The undeformed input preserves the actual failed
geometry. A geometry result must pass fresh final-SCF, spin, topology and
stationarity review, then receive its own new high-level electronic energy and
complete native Hessian/derivative/minimum/nuclear review. Old-geometry energies
or nuclear terms cannot transfer. The two [neutral signed branches](../phenothiazine_geometry_2026-11-02/README.md)
continue separately; no branch is selected before matched free energies exist.

Preflight exposed a portability error: the archived execution mount hid
current-checkout validation evidence referenced by the original neutral plans.
Those two attempts failed before Popen and have no scientific launch, SCF,
residency guard or remote writer. Their uncertain claim records and original
inputs remain intact. Exact old admission owners were confirmed closed; fresh
v2 namespaces capture byte-identical evidence outside the hidden checkout.
The [three namespace probes](private_namespace_preflight.json) verify every
dependency inside the actual archive mount without starting CUDA. The
[recovery record](visibility_recovery_admission.json) and
[live owner observation](controller_observation_01.json) keep preparation,
coordination and scientific starts separate. At that observation both controllers
were alive with exact recorded start times and one CPU each; all three scientific
cases were still unstarted and waiting for suitable capacity. No running
scientific worker was signaled or retried.
The later [actual start observation](actual_starts_01.json) and geometry reviews
above supersede those earlier waiting counts without changing their records.

The [116 focused checks](all_checks.xml) cover physical failure closure,
resource/residency interference, exact charge/spin and method profiles,
archive-visible byte capture, queue contracts, context-matched native sampling
and numerical conformer comparisons. They are software checks, not replacement
quantum validation. Use `python -m redox.validation.stationarity_failure --help`
to inspect the failed-gradient diagnosis and
`python -m redox.workflows.route_geometry --help` for preparation and execution.
The shared closed-failure reviewer also passes a new [actual neutral replay](neutral_shared_failure_review_01.json):
nine saved SCF frames, the same -25.47-cm-1 mode and maximum sampled-column error
7.75e-6. This rechecks saved evidence without starting a new quantum calculation.

The files `electronic_snapshot.json`, `native_snapshot.json` and
`sampling_snapshot_audit.json` were captured here while other queues continued.
They are **phthalimide polymer method controls, unrelated to phenothiazine**.
Their byte-identical canonical copies and scope are recorded in the
[corrected sampling publication receipt](../corrected_sampling_2026-11-02/snapshot_relocation_02.json).
Earlier mutable primary/neutral entry records are retained here as
`*_previous_current.json`; their owners and paths describe those earlier observations.
