# Session handoff — 2026-10-11, 20:08 Chicago results snapshot

This handoff supersedes stale progress descriptions in CONTINUATION_2026-10-11.md.
Read live records again before resuming; background jobs continue to change files.

## User objective and standing authorization

Work in parallel on physical validation of ALL computed properties and expanding
quinone/imide/etc motifs grafted onto the Merrifield monomer. Repair physical failures,
not merely label them. Make new-motif grafting/computation automatic and resumable.
Deliver verified results for completed candidates first; do not present preliminary
all-39 values as verified. Use all available suitable GPUs, without an arbitrary cap;
exclude occupied devices and cards too small for the actual workload. Update slides
as verified results arrive. Read AGENTS.md and the lambda-fleet skill before submissions.
Preserve unrelated existing dirty work; no project commit/push was made.

## Current deliverables

- Results: results/today_2026-10-11_verified_collection/LATEST_verified_candidates.xlsx
- Original proposals CSV: results/today_2026-10-11_verified_collection/LATEST_original39.csv
- Live feed with immutable snapshot pointers: results/today_2026-10-11_verified_collection/latest.json
- Checked 22-page slides: slides/redox_workflow.pdf
- Raw calculations: calcs/dft/<candidate>/<state>/
- Fleet queue, attempts, claims, logs: results/validation/candidate_fleet_today_2026-10-11/

Snapshot 20261011_200840_931573 reports 18/39 complete admitted native redox-potential
ladders, 114/115 original DFT geometry states, 84 admitted native states, 235 independent
arithmetic checks with zero discrepancies, and 7458 verified workbook cells.
This is numerical validation under the documented spin/solvent/thermal model, not
experimental accuracy. No exact graft has experimentally validated properties here.
Geometry completion, pipeline done markers, and complete property validation differ.

All 39 proposals plus one additional terminal-alcohol NDI are registered: 118 states,
351 checked conformers. Parent literature values are never assigned to grafted molecules.
Completed potential rows include methyl/ethyl/isopropyl NDI, AQ 1-O, methoxy BQ,
ether viologen, phthalimide, naphthalimide, dimethoxy phenothiazine, and nine D3TaLES grafts.
For example, ethyl NDI gives E1/E2 -0.7362/-1.4372 V vs Fc; methoxy BQ gives
-0.5702/-1.7853 V; ether viologen gives -0.6032/-1.1681 V. These are model predictions.

## Accuracy repairs achieved

Native solution-DFT Hessian/qRRHO/spin corrections replace incompatible gas-GFN
thermal transfer for admitted results. Conservative gradient finite differences,
full integration-grid response and matching CPHF grids validated singlet, radical,
and metal controls. Strict gates check source hashes, finite values, Hessian symmetry,
gradient, minimum frequencies, FD agreement, spin contamination and matching definitions.
Missing/failed/stale thermal records cannot silently fall back to legacy values.

Both Fc saddle geometries were repaired to native DFT minima. Fe occupation seeding
avoids an excited SMD solution about 0.944 eV too high. Current native Fc absolute
reference is 4.37121155 V. AQ-amino native ladder is -1.23177143/-2.14194378 V vs Fc.
d3_80izcq/red2's UMA-selected triplet was wrong: active DFT singlet is 1.7571 eV lower;
the correction, source archives and persistent spin override were adopted.
Other small/unknown even-state spin gaps receive explicit alternative DFT checks;
near-degenerate states need matched native spin thermals and ensemble treatment.
Large UMA gaps retain a preliminary spin-model assumption. Do not claim every
alternate spin or conformer has been independently exhausted.

Reorganization now has native SMD frozen-slow response controls and independent
CPU/GPU agreement. Corrected AQ-amino graft neutral and reduced endpoints both
needed genuine minimum repairs. Six fresh SCF legs give directional lambda
1.079353/1.191924 eV, mean 1.135639 eV, asymmetry 9.9126%; 44 independent checks pass.
See results/validation/physical_reorg_2026-10-11/graft_repaired_solution_pair/result.json.
This remains a qualified R2-geometry/WB97-energy composite with unequal-curvature
limitations, not a validated experimental rate. Exact 2,6-dimethyl-BQ control gives
1.391922 eV versus a measured fit interval 1.26–1.28 eV. Remaining flagged discovery
reorganization cases and conformer uncertainty still require work. The invalid
decomposed ion-pair helper stays excluded. Historical ~0.22 V benchmark error does
not establish accuracy of the revised native method; reassess matched benchmarks.

## Automation and fleet

Grafting/registration validates graph, stereo, charge and electron parity. Global
manifest writes serialize safely. Caches bind source geometry, spin, protocol and
checkpoint; repairs archive stale artifacts. Pipeline resumes preserved stages;
collector independently reconstructs arithmetic from raw sources and withholds
failed/pending properties. Latest completed full test run: 143 passed; rerun after
later edits, not as proof of physical accuracy.

Fleet expanded to 31 candidate assignments plus 2 physical-validation assignments
across lambda1/2/4/11/13. This is the launch allocation, not a live utilization count.
Recheck actual processes/reservations; finished/failed queue work must not be counted
as productive busy GPUs. 28 candidate done markers exist, but only 18 complete
potential ladders pass collection. Eleven unresolved failure markers currently exist:
d3_80prff, exp_ndi_n_butyl, d3_05ujsa, design_pmdi_methoxyethyl,
design_ndi_methoxyethyl, d3_05mnmt, design_ndi_methoxyethoxyethyl,
exp_viologen_methoxyethoxyethyl, exp_pmdi_trioxyethyl, exp_ndi_hexyl_pentyloxy,
d3_05lcla. Inspect each log for cause; do not treat preserved-partial refusal as
a second independent physical failure. Some old failure markers coexist with done
markers; done/current source evidence takes precedence.

Lambda13 H100 admission required an isolated matching NVIDIA565 userspace runtime;
system files/drivers were untouched. Use only on lambda13:
results/validation/lambda13_runtime_shim_2026-10-11/with_runtime.sh
and FAIRCHEM_CACHE_DIR=/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem.
CUDA, UMA, active SP and native thermochemistry replay controls passed. GPU1 there
has another user's Julia process and must remain excluded. V10032GB/H100 cards used;
16GB full-workflow cards excluded for measured memory limitations. All compute
workers have disjoint two-CPU affinities and atomic GPU reservations.

## Processes actually checked at handoff

- lambda1 PID686270: native physical cohort, GPU0/CPU0–1; still running.
- lambda1 PID935999: spin queue controller CPU4; actual child PID1940752 on GPU1,
  CPU2–3, design_pmdi_methoxyethoxyethyl:red2. Controller waits without reserving GPU.
- lambda1 PID2129640: bounded d3_80prff neutral repair, GPU7/CPU22–23, two-hour limit.
- lambda2 PID1296500: verified-results collector, CPU76, every120seconds.
- lambda2 PID1232655: checked-slide publisher, CPU77.

Collector/publisher stop at 23:55 Chicago. Candidate workers stop taking new work
at 23:30 Chicago; existing stages also have bounded deadlines. In a new session,
verify PID/start/command identity and source freshness before any signal/resume.
Do not launch duplicate controllers or bypass reservations.

## Active difficult repair: d3_80prff neutral

Actual full-response Hessian is finite but asymmetric by 0.004691949 Ha/Bohr²,
well beyond the unchanged 1e-6 gate; apparent minimum frequency about -17 cm^-1.
Tightening CPHF 1e-6 to 1e-10 did not change it. Isolated SMD non-electrostatic CDS
curvature reproduces the asymmetry: sharp SASA/contact behavior around H28.
Step reduction improves but does not pass the gate. Diagnostics are preserved.
Two signed H28 displacements followed by conservative R2-D4/SVP/SMD refinement
are being tested privately. The first reached120steps with force oscillations;
native/FD checks are in progress. Canonical geometry/SP remain unchanged and
failure marker must stay closed until a strictly validated replacement is adopted.
Any adoption needs SHA archives and recomputation of affected SP/cross points.
Do not blindly symmetrize, relax the gate, or change the solvent model to admit it.
PySCF2.14 geomopt.optimize(assert_convergence=True) asserts SCF convergence, not
GeomeTRIC optimization convergence; future repair reporting should capture
geomopt.kernel()[0] explicitly in addition to independent physical gates.

## Next session priorities

1. Read latest.json and inspect actual active workers; preserve running calculations.
2. Complete d3_80prff bounded repair and independently inspect other failure logs.
3. Continue core native minima and source-addressed alternate-spin/ensemble checks.
4. Advance remaining solution reorganization cases only with qualified current endpoints;
   assess conformer/cutoff/solvent/benchmark uncertainty before claiming reliability.
5. Resume admissible unfinished candidates using all suitable available GPUs through
   lambda-fleet; reuse verified caches, preserve failed history, never substitute parent values.
6. Refresh audited workbook and checked slides as new properties pass their gates.

Scientific Python: /nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python.
The Git worktree has substantial pre-existing unrelated modifications/deletions;
preserve them. Chat stream-disconnection cause was not established; saved results
and the background processes listed above were verified independently.
