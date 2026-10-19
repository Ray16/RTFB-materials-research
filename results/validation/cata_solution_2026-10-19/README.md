# Cata-imide solution-state calculations

Two new exact grafts, `hyp_cata_ndi_ethyl` and `hyp_cata_adi_ethyl`, extend the
imide search to adjacent imide annulation and a larger fused core. Their source
parents and CH2Cl2 measurements remain separate; no measured value is transferred.
See the [source reconstruction and UMA records](../../discovery/cata_core_expansion_2026-10-19/README.md).

Both designs are now registered in `config/custom_candidates.py` with six states
in the shared manifest. Canonical conformers are byte-identical to the isolated
build; all six [UMA warm starts](../../discovery/cata_core_expansion_2026-10-19/registration_prepared/UMA_adoption.json)
were checked and installed without recalculation. Original records remain intact.
The prior registry and byte-exact manifest are archived alongside that receipt.

The [plan](plan.json), [node preflight](lambda13_preflight.json),
[launch](launch.json) and [actual-start verification](start_verification.json)
separate preparation, resident predictors and completed records. Six independent
state workers use six suitable free H100 NVL GPUs, disjoint two-CPU sets, capped
native pools and atomic system reservations. The 24-GiB Hessian response budget
is independently qualified; the scientific protocol and physical gates are unchanged.

Each state runs full-response solution geometry optimization, the active
wB97M-V/TZVPD/SMD energy and derivative-checked native thermochemistry with
internal-mode projection. The canonical directories are
`calcs/dft/<candidate_id>/<neu|red1|red2>/`; per-state logs and ownership are under
`tasks/`. Each case has one actual attempt. Failed or partial outputs remain for
review; an observation timeout is not a reason to rerun it. Current producer/control
hashes are frozen for the running epoch.

The cata-NDI dianion's 0.2348-eV UMA gap requires a subsequent source-bound DFT
spin comparison. A completed selected-spin state alone cannot admit that second
wave. Near DFT spin gaps additionally require native spin thermals/ensemble review.
No complete solution potentials are claimed at launch. Geometry/SCF convergence,
stationarity, Hessian derivatives, minimum, integrity and spin admission remain
separate gates. Solution kinetics and exact-graft experimental performance remain
unvalidated.

The [prior-review epoch check](prior_review_live_epoch_check.json) proves that
new registration postdates the immutable full table review. The candidate feed
now contains 27 additional designs with pending quantities blank. Its current
export has 24 complete potential pairs and passes 576 arithmetic checks. The
[software suite](software_checks.log) passes 540 tests.

The preceding goal turn made progress: it published model-limited candidate
exports, curated four source graphs and completed six reserved UMA states,
pushed in commits `144cf51` and `1f8cfaa`. This epoch advances the two prepared
structures into the shared registry and the same solution protocol as prior cases.
