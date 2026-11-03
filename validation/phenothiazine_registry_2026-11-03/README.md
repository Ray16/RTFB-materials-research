# Qualified phenothiazine property in the shared candidate pool

The exact N-(4-methylbenzyl)phenothiazine local graft now has a
[verified candidate record](../../results/candidates/verified_cases/phenothiazine/current.json)
and an evaluation in the [current unified pool](../../results/candidates/pool/current.json).
Its first oxidation is +0.626086863 V versus Fc/Fc+ in MeCN at 298.15 K and the
50-cm-1 rotor cutoff. Both retained neutral branches, all three rotor settings,
state/reference components and their numerical policies remain available.
The numerical scenario bounds are not an experimental confidence interval.

The [registry comparison](comparison_02.json) proves that all 912 graphs, aliases,
declarations, route decisions and earlier evaluations are preserved. Only the
phenothiazine graph receives one appended evaluation. Graphs with reported
evaluations increase from 76 to 77; this count includes historical reports and
does not mean 77 fully qualified materials. Route counts remain three proceeds
and 909 holds, with zero full material admissions.

Grafting free energy, relaxation, stability beyond the declared states, coupling,
polymer occupancy and accessible-capacity evidence are not supplied by this
potential. Uncomputed property fields remain null. The all-axis decision stays
pending because the application criteria and remaining evidence are unresolved.
The eight substituted phenothiazine proposals inherit neither this value nor
the baseline route pass.

The maintained publisher is
`python scripts/pipeline/publish_qualified_redox_case.py --help`. It rechecks
the original redox assembly's source hashes and replays both exact cycles before
publishing. The candidate packet contains frozen redox, registry and policy
inputs. `redox.screening.qualified_cases.read_case` reconstructs it from those
local snapshots; it does not repeat raw quantum QC.

`config/candidate_registry_sources.json` now declares this qualified case as a
separate source. Normal registry rebuilding retains it automatically, using
the same case reader and exact-graph join. Frozen case inputs are copied into
the registry snapshot, so the join remains reviewable without private runtime
paths. A property cannot override a current route hold or create a material
admission. [Thirty-six checks](checks_02.xml) cover the new case, registry and
feasibility boundaries, including changed identity, values, method/reference,
numerical policy, component spin, model scope, missing fields and corrupt snapshots.

The manuscript already contains this value. Its frozen identity inventory and
scientific results are unchanged by this catalogue update.

A [relocated packet replay](relocated_case_replay.json) succeeds with original
candidate, manuscript, validation and raw-calculation paths blocked. This
replays the case arithmetic and field definitions; it does not repeat physical
quantum validation. The initial case prototype is retained on shared storage;
the context-complete v2 record is the maintained publication.
