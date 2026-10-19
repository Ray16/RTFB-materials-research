# Computed-state qualification before experimental calibration

Experimental source quality is distinct from computed-state quality. The reusable
state spin reader now lives in `redox.validation.spin_evidence`, with a compatibility
export for existing collectors. `reference_admission.event_admission` requires
current native endpoints, a matching live Fc/thermal definition, stated completed
spin admission, and agreement with the serialized potential. A live spin writer,
stale potential or unfinished close-spin decision cannot calibrate uncertainty.

The general benchmark scorer and scorecard uncertainty inputs now use this gate.
The independent audit reconstructs its own per-event decision from raw identities,
spin reports and partitions, including the same explicit UMA/standard-doublet
models. It does not import the primary admission predicate. Bad experimental
agreement is not an exclusion criterion: physically admitted outliers remain in
statistics. The +0.547 V DMeOEPT second-wave discrepancy is retained.

The full suite passes 529 tests. `table_review.log` records the fresh complete
cross-reader table verification; its status must pass before promotion. Source
changes were made only after all parent/imide spin workers were verified terminal.
