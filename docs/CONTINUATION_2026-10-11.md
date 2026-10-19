# Physical validation and candidate continuation — 2026-10-11

The property arithmetic and the physical calculations have separate completion gates.
The preserved `review_2026-10-11_spin_and_thermal_checkpoint` passed 10,223 independent
checks at its input snapshot. The reference and thermal records are now being repaired;
a fresh export must pass the independent audit after those producers finish.

The thermal implementation has completed independent conservative-gradient Hessian
controls for an organic singlet, organic radical and metal reference. Full integration-grid
response and matching CPHF grids resolved the analytic/numerical discrepancies. The
126-column AQ-dianion numerical reference has no imaginary modes; the original gas-GFN2
thermal correction differs from it by 0.377 eV. Although all 16 original imaginary-mode
states now have stable solution-model counterparts, these low-level minima alone cannot
validate transfer to the DFT free-energy surface.

Both original Fc geometries have genuine negative DFT modes: Fc+ at about −22 cm⁻¹
and neutral Fc at −35 cm⁻¹. Two signed-mode Fc+ refinements produce stable solution
minima with practically degenerate energies; the adopted state and its native thermal
record are validated. Neutral Fc refinement is running. The metric worker retains
complete archives and hashes. The
Fe occupation policy also now prevents an excited SMD solution approximately 0.944 eV
too high at the refined geometries. Matched native thermal records for both Fc states and
each candidate ladder are required before native-definition potentials can be exported.
Details are in the [thermal closure report](../results/validation/physical_thermal_spin_2026-10-10/THERMAL_PROTOCOL_CLOSURE.md).

The separate reorganization worker is running three controls followed by the exact
2,6-dimethylbenzoquinone experimental comparison. Its native SMD frozen-slow calculation
retains the actual SMD cavity, dielectric constants and CDS bookkeeping. Independent
CPU/GPU reference agreement is 1.07e-10 eV. The measured benchmark is 1.26–1.28 eV with
documented electrochemical-fit limitations; this is a useful matched comparison, not an
error-free standard. The 13 flagged discovery steps follow only after control qualification.
See the [accuracy contract](../results/validation/physical_reorg_2026-10-11/ACCURACY_CONTRACT.md)
and [launch manifest](../results/validation/physical_reorg_2026-10-11/native_smd_launch_manifest.json).

Candidate expansion contains 39 proposals representing 38 connectivities. Two pilots,
ether viologen and phthalimide, have five converged, graph-checked UMA warm starts and
no production solution DFT properties yet. AQ 1-O and terminal-alcohol NDI are now
registered too: their six states and 16 conformers pass an independent registration
audit, and all 202 previous manifest rows remain unchanged. The four registered motifs
have 11 declared electronic states. Their
[production queue](../results/discovery/candidate_campaign_2026-10-11/campaign.json)
uses the corrected thermal method. Further quinone/imide candidates remain prepared in
[follow-up inputs](../results/discovery/followup_inputs_2026-10-11/batch_config.py).
Synthetic witnesses and proposed resin-transfer steps remain separate from computed
properties. No literature parent metric is assigned to a grafted derivative.

The automation repair includes global manifest serialization, complete UMA cache
identities, explicit stereo preservation, route-note exports and spin-override cache checks.
The project test suite passes 125 tests at this code checkpoint. Strict readiness
now rejects thermal REVIEW states. Invalid or incompatible Fc thermal
records cannot define referenced potentials. New pipeline plans include the checked native
thermal calculation for both the motif and Fc before readiness. Run manifests distinguish
audited arithmetic, selected-state readiness and experimental accuracy; a passing audit
does not mark experimental accuracy validated. GPU workers remain capped at two, with
disjoint CPU allocations and shared fleet reservations.

## Expanded fleet and delivery

The user reaffirmed the standing instruction to use all available suitable GPUs.
`AGENTS.md` now records it. The original 39 proposals and the additional NDI are all
registered: 118 states and 351 independently checked conformers; manifest 315 rows.
The fleet now has 31 candidate workers plus two physical-validation workers.
It uses 26 V100 32 GB GPUs on lambda1/2/4/11 and seven H100 NVL GPUs on lambda13.
Occupied GPUs and the 16 GB cards unsuitable for the full workflow remain excluded.

Lambda13's driver-library mismatch was resolved for these jobs with an isolated
matching-userspace runtime; no system installation or driver reload occurred. The
occupied H100 remains excluded. A reserved numerical control verifies UMA to about
1 micro-eV, active SMD/gas single points to 1.3e-10 eV, and native thermochemistry
to 4.7e-8 eV across V100/H100. The shared model cache is explicitly addressed by its
NFS path. The [fleet manifest](../results/validation/candidate_fleet_today_2026-10-11/launch_manifest.json)
records assignments, dependencies and the private runtime evidence.

The repaired native Fc pair is complete and its occupation control passed; the live
reference is 4.37121155 V. The AQ-amino native ladder gives -1.23177143 and -2.14194378 V
vs Fc, with 0.91017235 eV disproportionation energy. These are validated calculations
under the documented model, not measurements of the grafted structure. The 49-state
core thermal campaign continues, including true minimum repairs and dependent cross points.

New-candidate small or unknown even-state spin gaps require additional DFT checks.
Near-spin states remain withheld pending matching native thermal/ensemble treatment.
A source-addressed spin controller waits without occupying a GPU until suitable
candidate electronic records and a free validation slot are available.

The [verified-results feed](../results/today_2026-10-11_verified_collection/latest.json)
updates immutable CSV/Excel snapshots every two minutes through 23:55 Chicago.
It distinguishes all original 39 proposals from the additional NDI and independently
checks raw inputs/arithmetic, current spin evidence and changing-source hashes.
Pending or flagged primary values stay blank; model-surface estimates are separate.
The slide publisher watches this feed, generates result frames, builds and checks the
PDF, and publishes `slides/redox_workflow.pdf` after completed candidate potentials arrive.
The existing methods were updated and earlier figures moved to a labeled historical appendix.
Original slide files are preserved under the fleet validation directory.

This checkpoint supersedes running-job descriptions and the earlier voluntary two-worker
limit in the previous progress report. Field-specific completion must be read from the live
feed and actual workers; no full new-candidate result is claimed from geometry completion.
