# Candidate results

**Corrected thermal protocol admitted (2026-10-19).** Current exports require
explicit internal-mode thermochemistry and a consistently corrected Fc reference.
The implementation is integrated into the producer and both primary readers;
232 cases passed 6,961 independent qualification checks. Genuine soft internal
modes are retained. Earlier unprojected snapshots are historical predictions.
See the [correction and evidence](../validation/rigid_mode_thermal_2026-10-18/README.md).
Arithmetic agreement establishes reproducibility, while each reported field still
requires its own geometry, derivative, spin and integrity QC. Solution kinetics,
durability and exact-graft experimental accuracy remain unvalidated.

Open [latest.xlsx](latest.xlsx) and select the **Additional** sheet for the new
candidates. [latest.csv](latest.csv) contains just those additional candidates,
and [comparisons.csv](comparisons.csv) contains shifts from matched exact-graft
controls. These stable paths update with the active campaign exports.

The active batch contains 25 additional candidates. The workbook and current
snapshot report which potential pairs are complete. Blank numeric fields are pending or excluded by
field-specific QC. These are predictions for exact molecular structures; parent
measurements are separate and exact-graft experimental accuracy is unvalidated.

The immutable snapshot and its evidence are indexed by [current.json](current.json).
Exact structures and candidate IDs are included in the CSV/workbook; definitions
are in [config/custom_candidates.py](../../config/custom_candidates.py). Raw
geometry, energy and native-thermal records live in
`calcs/dft/<candidate_id>/<state>/`, with the initial UMA records in
`calcs/uma/<candidate_id>/<state>/`.

The checked slide deck is [slides/redox_workflow.pdf](../../slides/redox_workflow.pdf).

Selected spin equilibria are explicit in the state sheet: the selected-surface
native correction, effective correction used for G, equilibrium mixing term and
alternative population are separate fields. Single-surface kinetic validation is
withheld for mixed-spin states. Integrity must pass for every included component
pair. The current batch has 24 complete potential pairs out of 25; consult the
feed timestamp for any later updates.
