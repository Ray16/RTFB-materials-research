# Maintained polymer figures

| Figure | Evidence and scope | Plotting script |
| --- | --- | --- |
| [Accepted spacing conformers](spacing_conformers.pdf) | Actual phthalimide XYZ snapshots, core contacts and within-graph G differences; four accepted local minima in finite two-pose sets | `scripts/plotting/polymer_models/plot_spacing_conformers.py` |
| [Finite ensemble cycles](conditional_ensemble_baseline.pdf) | Twelve native states, two-pose nonadditivity and prescribed rotor-cutoff sensitivity | `scripts/plotting/polymer_models/plot_conditional_ensemble_baseline.py` |
| [Exact spacing connectivity](spacing_construction.pdf) | Open five-repeat phthalimide graphs, styrene padding and residual benzyl chlorides | `scripts/plotting/polymer_models/plot_spacing_construction.py` |
| [Exact PMDI construction](pmdi_occupancy_construction_compact.pdf) | Candidate occupancy graphs; preparation alone supplies no energy or synthesis eligibility | `scripts/plotting/polymer_models/plot_candidate_occupancy_construction.py` |
| [PMDI capacity ceilings](pmdi_capacity_bounds.pdf) | Composition/scaffold bounds with stated electron budgets and utilization assumptions | `scripts/plotting/polymer_models/plot_polymer_capacity_bounds.py` |
| [Segment deformation terms](segment_deformation_components.pdf) | Signed sampled-fragment reference differences, including cancellation | `scripts/plotting/polymer_models/segment_deformation_components.py` |

The [spacing-conformer data](spacing_conformers.json) include graph-ordered
coordinates, core groups, explicit molecular identities, proper view rotations,
native G parents and source hashes. The [cached replay](spacing_conformers_replay.json)
has identical pixels while an audit hook blocks original calculation, review
and campaign inputs. This reproduces the drawing; it does not repeat physical
qualification. Reproduce in the plotting environment with a fresh output path:

```bash
PYTHONPATH=src python scripts/plotting/polymer_models/plot_spacing_conformers.py \
  --cached results/figures/polymer_models/spacing_conformers.json \
  --out /tmp/spacing_conformers_replay/spacing_conformers
```

Apply the documented one-thread native pool caps. Original-source regeneration
uses `--review results/reviews/2026-10-23/baseline_native_complete/summary.json`.
The manuscript's frozen figure copies preserve their own evidence epoch.

The [local motif contact plot](local_motif_contacts.pdf) shows the four checked phthalimide capped-fragment electronic interactions versus sampled core separation. Cap-length ranges are not uncertainty intervals; basis-superposition correction and conformer convergence remain open. See the [source and checks](../../reviews/2026-10-26/local_motif_contact_figure/README.md).
