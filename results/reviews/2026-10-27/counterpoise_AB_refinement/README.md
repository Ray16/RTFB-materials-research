# Full-complex counterpoise checkpoint refinement

The five gas electronic terms finished their production calculations. Independent
checks accepted A-own, A-ghost, B-own and B-ghost. The original AB checkpoint
reproduced its energy exactly, but its orbital-gradient norm was 1.13037e-7,
above the unchanged 1e-7 acceptance threshold. That attempt remains unqualified
and preserved.

A separate eight-thread attempt starts from that exact reviewed AB checkpoint,
with SCF stopping targets tightened to 1e-13 hartree and 1e-9 for the orbital
gradient. Geometry, common quadrature, method, auxiliary basis and independent
acceptance threshold are retained. The [startup record](startup_verification.json)
checks the actual process identity, CPU binding and thread pools.
The [independent-review observer](observer_launch.json) waits for completion,
then checks the checkpoint and recomputes the orbital gradient. Its live status
is `calcs/validation/2026-10-27_AB_polish_review/status.json`.

The full gas counterpoise interaction remains unqualified until all five terms
have matching independent stationarity evidence. This does not supply a MeCN
correction, association free energy or polymer ensemble. Raw refinement outputs
are in `calcs/local_motif_counterpoise/cpu_one_spacer_AB_polish_2026-10-27/`.

The [explicit collection manifest](collection_manifest.json) joins the new AB attempt with the four already checked fragment terms; its [finite collector](collection_launch.json) waits for the independent AB check and then reuses the original five-term admission gates. Its live result is `calcs/validation/2026-10-27_full_counterpoise_checked_collection/status.json`. Preparing this handoff does not supply a counterpoise energy.
