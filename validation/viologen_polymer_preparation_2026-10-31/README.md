# Exact viologen: charged polymer preparation

The route-bound baseline viologen now has **12 audited initial polymer models**:
six occupancy graphs in each of two declared backbone configurations. The
[actual connectivity figure](../../results/figures/polymer_models/viologen_occupancy_construction.pdf)
compares adjacent sites 2/3 with sites 2/4 separated by an unreacted benzyl-chloride
site. All models preserve the open five-repeat chain and plain-styrene padding
at sites 1 and 5. This preparation advances the exact-candidate polymer setup;
it does not complete reviewer M1 or the occupied-ensemble requirement M4.

The [independent replay](checked_replay.json) verifies selected XYZ connectivity,
stereochemistry, covalent bond lengths, benzyl--pyridinium N attachment,
remaining covalent chlorides, and four atom/charge-balanced cycles. Every
backbone has ten carbons, nine backbone C--C bonds and two ends. Analytic
quaternization masses agree with masses of the actual graphs plus the declared
counterions. A graft releases chloride without removing precursor hydrogen.

The [preparation manifest](../../results/polymer_models/viologen_preparation_2026-10-31_v2/preparation.json)
and adjacent `plan.json`, source snapshot, per-case metadata, trial XYZs and
embedding diagnostics preserve inputs and partial progress. The first attempt
stopped after four graphs when ETKDG yielded no embedding; its original outputs
remain under `results/polymer_models/viologen_preparation_2026-10-31/`. The second
attempt adds one bounded random-coordinate fallback and visible per-case failure
records. Four graphs need that fallback; another yields only one of four requested
starts. Forty-five trial geometries are retained. These are not a complete
conformer search. The 200-step MMFF pre-relaxations did not converge; their
energies rank starting coordinates only and never enter solution thermodynamics.

Coordinates contain the bare polymer, with declared omitted separated spectator
ions. The fully oxidized **PF6 dry-salt mass reference is separate**: two PF6
anions per viologen. One/two grafts have reference masses 1093.886/1519.580 g/mol
and assumed two-electron capacity ceilings 49.00/70.55 mAh/g. Neither ceiling
establishes accessibility, achievable loading or mobile-ion mass during redox.
Precursor iodide, released chloride and later PF6 exchange are distinct species
and stages. Association/exchange free energies have not been computed.

No UMA, xTB, SCF, native Hessian or GPU job was started for this preparation.
The manifest keeps charged-solver/resource admission, solution geometry,
conformer sampling, native G and material admission false. Global charge/spin
hypotheses do not establish spin ordering or localization. The existing neutral
calculation entry point still rejects charged graphs. Current live neutral
polymer calculations are separate; new GPU starts remain subject to the shared
storage-headroom gate.

## Reproduce

Use the redox environment with native CPU thread pools capped at one. The
preparation/check commands allocate one disjoint CPU and preserve prior outputs.

```bash
python scripts/pipeline/prepare_charged_polymer.py \
  --config config/viologen_polymer_preparation.json \
  --registry results/candidates/pool/2026-10-30_viologen_checked/registry.json \
  --routes results/candidates/pool/2026-10-30_viologen_checked/inputs/route_evidence.json \
  --assessment results/candidates/pool/2026-10-30_viologen_checked/inputs/assessment.json \
  --out /path/to/new/preparation

python scripts/validation/check_charged_polymer_preparation.py \
  --preparation results/polymer_models/viologen_preparation_2026-10-31_v2/preparation.json \
  --out /path/to/new/replay.json

python scripts/plotting/polymer_models/plot_viologen_occupancy.py \
  --preparation results/polymer_models/viologen_preparation_2026-10-31_v2/preparation.json \
  --audit validation/viologen_polymer_preparation_2026-10-31/checked_replay.json \
  --out /path/to/new/figure \
  --renderer-python /path/to/python/with/PyMuPDF
```

The construction/audit suite passes 53 tests, including actual-geometry fallback,
total-embedding failure, failure persistence, exact graph masses, ion/spin
rejection, changed capacity, unbalanced cycles and forbidden qualification.
The new method and figure are included in SI Section 3.4, with the preceding
local-grafting and molecular results preserved. The paper remains a working draft.
