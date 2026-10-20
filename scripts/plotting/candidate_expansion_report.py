#!/usr/bin/env python
"""Render one data-linked presentation of design coverage and existing evaluation."""

import argparse
import hashlib
import html
import json
from pathlib import Path

from redox.core.io import atomic_write_json


def matched_figure(out, review):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    cores = ("viologen", "pmdi", "ndi")
    axes = ("cap_length", "linker_C2", "ring_F", "ring_CN", "ring_OMe", "ring_Me")
    labels = ("Cap length", "Ether linker", "Ring F", "Ring CN", "Ring OMe", "Ring Me")
    values = np.full((2, len(cores), len(axes)), np.nan)
    seen = set()
    for row in review["comparisons"]:
        if row["core"] not in cores or row["contrast"] not in axes:
            continue
        i, j = cores.index(row["core"]), axes.index(row["contrast"])
        if (i, j) in seen:
            raise ValueError(
                "Ambiguous family/axis figure cell; define a separate panel"
            )
        seen.add((i, j))
        for step, field in enumerate(("delta_E1_V", "delta_E2_V")):
            if row[field] is not None:
                values[step, i, j] = row[field]
    vmax = (
        max(float(np.nanmax(np.abs(values))), 0.01)
        if np.isfinite(values).any()
        else 0.01
    )
    fig, panels = plt.subplots(1, 2, figsize=(13, 4.1), layout="constrained")
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#eeeeee")
    for step, panel in enumerate(panels):
        grid = values[step]
        drawn = panel.imshow(grid, vmin=-vmax, vmax=vmax, cmap=cmap, aspect="auto")
        panel.set_xticks(range(len(axes)), labels, rotation=30, ha="right")
        panel.set_yticks(range(len(cores)), ["Viologen", "PMDI", "NDI"])
        panel.set_title(f"Step {step + 1}: candidate − declared control")
        for i in range(len(cores)):
            for j in range(len(axes)):
                value = grid[i, j]
                panel.text(
                    j,
                    i,
                    "pending" if not np.isfinite(value) else f"{value:+.3f}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white"
                    if np.isfinite(value) and abs(value) > 0.65 * vmax
                    else "#10293b",
                )
    fig.colorbar(drawn, ax=list(panels), shrink=0.85, label="Potential shift (V)")
    fig.suptitle(
        "Matched design effects • computed MeCN potentials\nSign and magnitude do not establish battery improvement",
        fontsize=13,
    )
    for extension in ("png", "pdf", "svg"):
        fig.savefig(
            out / f"matched_design_effects.{extension}", dpi=220, facecolor="white"
        )
    plt.close(fig)


def cell(value):
    if value is None:
        return "pending"
    if isinstance(value, float):
        return f"{value:+.3f}"
    return html.escape(str(value))


def table(headers, rows):
    return (
        '<div class="scroll"><table><thead><tr>'
        + "".join("<th>" + cell(x) + "</th>" for x in headers)
        + "</tr></thead><tbody>"
        + "".join(
            "<tr>" + "".join("<td>" + cell(x) + "</td>" for x in row) + "</tr>"
            for row in rows
        )
        + "</tbody></table></div>"
    )


def render(out):
    source = out / "current.json"
    review = json.loads(source.read_text())
    count = review["counts"]
    matched_figure(out, review)
    page = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Candidate design and evaluation</title><style>
body{font:16px/1.5 system-ui,sans-serif;color:#16314a;max-width:1500px;margin:30px auto;padding:0 20px}
a{color:#0072b2}.note{padding:16px;background:#fff4df}.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:14px}td,th{padding:9px;text-align:left;border-bottom:1px solid #ddd}th{background:#edf3f7}input{padding:10px;width:90%;max-width:600px}.flow{padding:16px;background:#edf3f7}td{vertical-align:top}</style>
<h1>Candidate expansion and evaluation</h1>
<p><a href="../../../docs/CANDIDATE_EXPANSION_LOOP.md">Workflow and reproduction</a> · <a href="../README.md">Discovery guide</a> · <a href="current.json">Exact data and provenance</a> · <a href="designs.csv">Design ledger</a> · <a href="comparisons.csv">Matched effects</a></p>
<p class="flow">Define a design question → select matched structures and controls → review exact grafting routes → evaluate with the existing pipeline → compare admitted fields → define the next question.</p>"""
    page += f"<p><b>{count['unique_designs']} distinct graphs after joining {count['enumerated_pool']} enumerated designs with {count['selected_candidates']} current candidates and {count['external_controls']} external control.</b> These sets overlap; no exact graph is counted twice. {count['matched_comparisons']} declared comparisons ({count['complete_matched_pairs']} with both potential shifts).</p>"
    page += (
        "<p>Candidate evaluation snapshot: "
        + cell(review["source_candidate_as_of"])
        + ". This collection joins existing results and submits no calculations.</p>"
    )
    page += f"<p class='note'><b>{count['prospective_route_proceed']} designs currently pass the prospective full-redox route gate.</b> Missing candidate-specific route evidence stays explicit. Historical computations remain available. This is a separate dimension from redox computation and does not exclude base-required BQ/AQ designs on their HCl-reference energy.</p>"
    page += "<h2>Coverage by declared core or source family</h2><p>Counts refer to unique exact grafts. Source-family labels are retained where a finer core assignment was not declared. These are coverage counts, not efficacy rankings.</p>"
    page += table(
        [
            "Core / source family",
            "Unique designs",
            "In enumeration",
            "Selected candidates",
            "External controls",
            "Potential pairs",
            "Route proceed",
        ],
        [
            [
                r[k]
                for k in (
                    "family",
                    "unique_designs",
                    "in_enumerated_pool",
                    "selected_candidates",
                    "external_controls",
                    "admitted_potential_pairs",
                    "prospective_route_proceed",
                )
            ]
            for r in review["family_coverage"]
        ],
    )
    page += "<h2>What changes relative to a declared control?</h2><p>Δ = candidate minus control. Both endpoints must supply the same charge/electron event under the existing evaluator. Pending cells remain blank in CSV. A larger voltage shift, narrower wave spacing or lower mass is not an automatic selection benefit.</p>"
    page += '<figure><img src="matched_design_effects.png" style="width:100%" alt="Signed first and second potential shifts for eighteen matched design contrasts"><figcaption>Common six-axis panels use the same declared family controls. <a href="matched_design_effects.pdf">PDF</a> · <a href="matched_design_effects.svg">SVG</a>. Additional source-guided comparisons are retained below.</figcaption></figure>'
    page += table(
        [
            "Candidate",
            "Control",
            "Design axis",
            "ΔE1 (V)",
            "ΔE2 (V)",
            "Δ spacing (V)",
            "Δ mass (g/mol)",
            "Field status",
        ],
        [
            [
                r[k]
                for k in (
                    "id",
                    "baseline",
                    "contrast",
                    "delta_E1_V",
                    "delta_E2_V",
                    "delta_wave_spacing_V",
                    "delta_MW",
                    "status",
                )
            ]
            for r in review["comparisons"]
        ],
    )
    page += "<p><a href='design_cells.csv'>Full core × anchor × cap × ring coverage</a> preserves unselected cells. All symmetry-equivalent generation recipes remain in the design ledger.</p>"
    page += "<h2>Selected candidates and evaluation feedback</h2><p>Polymer counts are checked electronic structures for the same attachment fragment, in different polymer models. They do not transfer polymer energies into the capped-site property columns. Model limitations and route admission remain separate.</p><label>Filter candidates <input id='filter' placeholder='ID, family, design axis or status'></label><div id='candidates'>"
    selected = [r for r in review["records"] if r["evaluation_role"] is not None]
    page += table(
        [
            "Candidate",
            "Family",
            "Axis",
            "E1 / E2 (V vs Fc)",
            "Redox stage",
            "Model validity",
            "Route / label",
            "Checked polymer structures",
            "Next evaluation question",
        ],
        [
            [
                r["id"],
                r["core"],
                r["design_axis"],
                f"{cell(r['E1_vs_Fc_V'])} / {cell(r['E2_vs_Fc_V'])}",
                r["evaluation_stage"],
                r["model_validation_status"],
                r["selection_label"] or r["route_status"],
                r["checked_polymer_models"],
                r["evaluation_feedback"],
            ]
            for r in selected
        ],
    )
    page += "</div><h2>Planned extensions</h2>" + table(
        ["Axis", "Status", "Question"],
        [[r[k] for k in ("axis", "status", "question")] for r in review["future_axes"]],
    )
    page += (
        "<h2>Interpretation limits</h2><ul>"
        + "".join("<li>" + cell(x) + "</li>" for x in review["limits"])
        + "</ul>"
    )
    page += """<script>document.getElementById('filter').addEventListener('input',e=>{const q=e.target.value.toLowerCase();document.querySelectorAll('#candidates tbody tr').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));});</script></html>"""
    (out / "index.html").write_text(page)
    atomic_write_json(
        out / "report_manifest.json",
        dict(
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            output_sha256=hashlib.sha256((out / "index.html").read_bytes()).hexdigest(),
            figure_sha256={
                path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in out.glob("matched_design_effects.*")
            },
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    render(parser.parse_args().out.resolve())
