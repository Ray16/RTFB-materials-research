# Integrated candidate expansion and evaluation

Open the [report](index.html) or [matched-effects figure](matched_design_effects.pdf).
The [workflow guide](../../../docs/CANDIDATE_EXPANSION_LOOP.md) and
[design plan](../../../config/candidate_expansion_plan.json) define how the existing
design space connects to the evaluator and how results inform the next batch.

Use [status](status.json) for counts and the source feed's date. The
[design ledger](designs.csv), [paired effects](comparisons.csv),
[family coverage](family_coverage.csv) and [full design cells](design_cells.csv)
have stable filenames. [Current JSON](current.json) binds these to the exact
input hashes, producer hashes, source evaluator snapshot and shared-storage
publication snapshot. Earlier snapshots remain in `snapshots/` on NFS.

This collection submits no calculations and changes no molecular registry.
Candidate-specific route review is a separate prospective gate. Historical
computed fields, scoped model limits and pending quantities remain visible;
base-required BQ/AQ designs are retained for evaluation.
