# Maintained candidate comparisons

Open the [comparison browser](index.html) for all candidate designs, sortable
properties, exact structures, filters and comparisons of up to six structures.
Open the local HTML from this collection folder; the interactive plot and table
data are embedded. All candidates remain in the table, including unfinished
calculations. Structure cards also show their proposed vinyl counterparts.

The collection uses five folders, with this guide and the browser at the top level:

| Folder | Contents |
| --- | --- |
| [plots/](plots/) | Three focused comparison plots, in PNG, PDF and SVG |
| [structures/](structures/) | Molecular plates and oxidation-state drawings |
| [supplementary/](supplementary/) | Complete candidate ladders and family coverage |
| [data/](data/) | Candidate tables, structure records and selection details |
| [metadata/](metadata/) | Captions, checksums, provenance and generated slide include |

| Stable figure | What it compares |
| --- | --- |
| [Potential–capacity map](plots/redox_capacity.png) · [PDF](plots/redox_capacity.pdf) | Admitted first/second potentials and nominal molecular capacities |
| [Synthesis tradeoff](plots/synthesis_tradeoff.png) · [PDF](plots/synthesis_tradeoff.pdf) | SA and nominal capacity, highlighting selected descriptor leaders |
| [New candidate shifts](plots/new_candidate_shifts.png) · [PDF](plots/new_candidate_shifts.pdf) | Matched contrasts with lower potentials for both events, one per core |
| [Molecular structures](structures/candidate_structures.svg) · [PDF](structures/candidate_structures.pdf) | Computed capped-model leaders and uncomputed vinyl-monomer pilot; attachment and vinyl group marked |
| [Redox-state schemes](structures/redox_states.pdf) · [computed-model SVG](structures/redox_states_computed.svg) · [vinyl-proposal SVG](structures/redox_states_vinyl.svg) | Three redox states per molecule; charge, multiplicity and electron-transfer steps |
| [All-candidate ladders](supplementary/redox_ladders.png) · [PDF](supplementary/redox_ladders.pdf) | Every candidate, with pending endpoints blank and nominal capacity alongside |
| [Family coverage](supplementary/coverage.png) · [PDF](supplementary/coverage.pdf) | Admitted potential calculations and pending/excluded cases by family |

The first three are the focused main views. The complete ladders and coverage are
supplementary views. [Highlights and selection rules](data/highlights.json) distinguish
the partial potential–capacity front from the existing qualified primary Pareto
policy. SA is an annotation; no missing kinetics or model uncertainty is filled
to claim a complete front. [Figure captions](metadata/captions.md) state methods and scope.

The modeled structures are p-methylbenzyl-capped local graft models, not the
vinyl-bearing monomers required for monomer-first polymerization. Computed
descriptor leaders therefore do not establish polymer synthesis priorities.

The redox-state companion uses illustrative formal resonance structures, checked
for valence, molecular charge, unpaired-electron count and preservation of atoms,
hydrogens and connectivity. Charges/radicals can be delocalized; these drawings
are not calculated density maps or measured bond orders. Capped-model potentials
are shown only where already admitted. Vinyl-state potentials and ground-spin
assignments remain uncomputed; multiplicities there are proposals. Arrows denote
the electron-count sequence, not demonstrated cycling reversibility.

Download the complete [CSV](data/all_candidates.csv) or [JSON](data/all_candidates.json).
The [manifest](metadata/manifest.json) records the audited input snapshot, candidate IDs,
coverage and output checksums. These filenames are overwritten on each successful
update. Earlier scientific snapshots stay in the existing immutable export archive.

The focused potential plot contains eligible two-electron anolytes only, with
three numbered candidates in consistent colors and other points in gray.
Catholytes remain in the complete table and supplementary comparison. Legends
are checked for collisions with axis text before export.

Create the isolated figure environment from the project's `redox` environment
(already prepared in this workspace). Visualization dependencies stay separate
from the live scientific workers:

```bash
python -m venv --system-site-packages .venv-figures
.venv-figures/bin/python -m pip install -r requirements-visualization.txt
```

To update from the current verified feed:

```bash
.venv-figures/bin/python scripts/plotting/candidates/update_comparison_collection.py
```

To collect newly finished results and update this collection plus the existing
historical comparison figures in one command:

```bash
.venv-figures/bin/python scripts/pipeline/collect_verified_candidates.py \
  --out results/reviews/2026-10-18/verified_candidates \
  --queue results/validation/aryl_solution_2026-10-19 \
  --update-figures
```

Add `--watch` to refresh after each audited collection while the queue runs.
These commands read existing scientific calculations; they submit no new jobs.
Plotly supplies the interactive plot and CairoSVG exports vector molecular
structures. The existing slide publisher's `--update-plots` option updates
this same collection through `update_candidate_overlays.py`.

Potentials use the current corrected native protocol and matched Fc reference.
Nominal capacity excludes solvent, supporting electrolyte, counterions and polymer.
Model window capacity is separately gated and remains provisional. Reorganization
fields stay blank when solution-level qualification is missing; gas descriptors
are not substituted. Source-parent measurements are separate from exact grafts.
These comparisons establish neither measured solubility and durability nor
experimental superiority. An admitted calculation is not complete scientific validation.

If scientific inputs change during rendering, the previous completed collection
is preserved and the command asks for a fresh verified snapshot. It does not
relax admission to fill a chart.

For monomer-first polymerization, a chemical route assessment should accompany
the plots. Conventional diimide and viologen branches are reasonable initial
families; exact asymmetric precursor preparation and grafting remain to establish.
Viologen reaction with 4-vinylbenzyl chloride has
[published precedent](https://pubs.acs.org/doi/10.1021/acspolymersau.1c00036), and
[styrene-functionalized NDI polymerization](https://www.mdpi.com/1422-0067/26/10/4807)
supports the broader scaffold concept, not these exact structures or battery use.
Quinones can [retard or inhibit radical polymerization](https://doi.org/10.1016/0014-3057(73)90026-8).
The BQ models' high nominal capacities therefore do not make them automatic
monomer synthesis priorities. These are chemical assessments, not verified routes.
