# Project highlights report

Rebuild from the pinned local evidence and PDF figure snapshot:

```sh
python scripts/report/build_report.py
```

Requires `dvisvgm` on PATH. Output: `results/report/project_highlights.html`.
The HTML is self-contained for viewing, including eight isolated SVG images;
source-document links and the hash manifest require the project workspace.
No external fonts, network requests or JavaScript packages are required to view it.
The report is a summary, not new physical validation. Update narrative and pinned
sources together. Missing figures fail the build rather than disappearing.

Browser validation (optional dependency: Playwright and Chromium):

```sh
python scripts/report/check_report.py
```

Checks cover all eight image decodes, unique IDs, anchor targets, numeric chart
switching, figure dialog/Escape, JavaScript errors and page overflow at 390,
768 and 1440 pixels. Screenshots and checks are written alongside the report.
Desktop hero and neighbor-chart screenshots were visually inspected during this
revision. Prior generator and HTML are preserved in `results/report/prior/`.

Source hashes live in `results/report/manifest.json`; its Git HEAD identifies the
build's starting repository state, while explicit generator/template hashes
identify the report source exactly. Bulk screenshots and HTML stay in runtime
storage; the generator, template, checks and evidence manifest are versioned.

## Visual shortlist revision

Install Python dependencies with `pip install -r scripts/report/requirements.txt`
and Chromium with `playwright install chromium`. The local validation environment
is `scripts/report/.check`; its browser cache is `scripts/report/.browsers`.

The generator now uses RDKit to depict three exact registered local-site graphs,
with the benzyl attachment bond selected by a unique substructure match. It embeds
an interactive archived Pareto comparison and a five-repeat PMDI composition
explorer. `shortlist.json` records evidence-based follow-up priorities, not a
validated best-material ranking. The current primary front remains empty.

`data/historical_pareto.csv` and `data/historical_scorecard.csv` are unmodified
snapshots from `/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/results/`.
The archived front flags are displayed as recorded, not recomputed or promoted.
Partial kinetics and thermal flags and the naphthoquinone route concern are
retained. Historical salt-mass capacities and the neutral PMDI fragment bound
are explicitly different accounting conventions.

Browser checks also cover molecular-image decoding, salt scenario switching,
and zero/full loading and one/two-electron capacity values. Desktop shortlist,
Pareto and capacity screenshots were visually reviewed. Source hashes capture
all input datasets and presentation code.
