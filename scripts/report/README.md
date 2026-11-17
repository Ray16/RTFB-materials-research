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
