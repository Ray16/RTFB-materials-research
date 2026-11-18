#!/usr/bin/env python
"""Figure 2 (manuscript, JCIM): the nine-stage screening pipeline.

A corrected, manuscript-accurate successor to
scripts/plotting/pipeline/plot_pipeline.py, which was stale for this figure
(only six boxes; no thermochemistry/properties/QC/scoring stages; and an
out-of-date "MAE 0.13 V (n=7)" validation gate that contradicts the maintained
benchmark). Validation is its own figure (Fig 4), so it is NOT drawn here.

This renders the actual spine of `redox.workflows.motifs.commands()`:
  build -> UMA pre-opt -> DFT geometry+SMD(MeCN) per state [neu/ox1/ox2]
  -> refined single-points -> thermochemistry (two-tier) -> derived properties
  -> QC coverage/cross-points -> scoring/Pareto,
with each stage tagged GPU/DFT vs CPU, and the level-of-theory / standard-state
conventions as an inset table (manuscript Table 1).

Needs the system `dot` executable (Graphviz).

  python scripts/plotting/manuscript/fig2_pipeline.py
    -> manuscript/figures/fig2_pipeline.png (+ .caption.md)

Physics, not fitting: this only draws the workflow; it changes no number.
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

import graphviz

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = ROOT / "manuscript" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

INK = "#111111"
MUTE = "#555555"
# compute-resource palette
FILL_GPU = "#CBD9EC"   # DFT / GPU stages      (blue family)
EDGE_GPU = "#2166AC"
FILL_ML = "#FBE3DE"    # ML (UMA) GPU stage    (warm)
EDGE_ML = "#D6604D"
FILL_CPU = "#E7EDE7"   # CPU stages            (muted green-grey)
EDGE_CPU = "#4D4D4D"
FILL_GATE = "#DCEBDD"  # QC / scoring accent
EDGE_GATE = "#1B7837"


def _stage(g, name, n, title, detail, fill, edge, tag):
    """A numbered stage box: index badge + bold title + detail + GPU/CPU tag.

    Every font is >= 18 pt (manuscript style contract: 18 pt floor everywhere)."""
    label = (
        f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="4" CELLSPACING="0">'
        f'<TR>'
        f'<TD VALIGN="MIDDLE"><FONT POINT-SIZE="26"><B>{n}</B></FONT></TD>'
        f'<TD VALIGN="MIDDLE"><FONT POINT-SIZE="20"><B>{title}</B></FONT>'
        f'<BR/><FONT POINT-SIZE="18" COLOR="{MUTE}">{detail}</FONT></TD>'
        f'<TD VALIGN="MIDDLE"><FONT POINT-SIZE="18"><B>{tag}</B></FONT></TD>'
        f'</TR></TABLE>>'
    )
    g.node(name, label=label, fillcolor=fill, color=edge)


def build():
    g = graphviz.Digraph("redox_pipeline_fig2", format="png")
    g.attr(rankdir="TB", bgcolor="white", nodesep="0.30", ranksep="0.42",
           splines="true", pad="0.3")
    g.attr("node", shape="box", style="rounded,filled,setlinewidth(2)",
           fontname="DejaVu Sans", fontcolor=INK, penwidth="2.4",
           margin="0.22,0.14")
    g.attr("edge", color="#3A3A3A", penwidth="2.2", arrowsize="0.9",
           fontname="DejaVu Sans", fontsize="18", fontcolor="#333333")

    g.node("title",
           label=('<<FONT POINT-SIZE="26"><B>Screening pipeline: motif '
                  '&#8594; reliability-gated redox properties</B></FONT><BR/>'
                  '<FONT POINT-SIZE="18" COLOR="' + MUTE + '">'
                  'one reproducible driver, nine content-addressed stages; '
                  'E&#176; referenced to Fc/Fc&#8314; in SMD(MeCN)</FONT>>'),
           shape="plaintext", style="", margin="0")

    GPU = "GPU / DFT"
    ML = "GPU / ML"
    CPU = "CPU"
    stages = [
        ("build", 1, "Build grafted-motif model",
         "RDKit ETKDG ensemble &#183; 4-methylbenzyl (Merrifield) cap<BR/>"
         "redox states per motif (q, M)", FILL_CPU, EDGE_CPU, CPU),
        ("uma", 2, "UMA MLIP pre-optimization (gas phase)",
         "charge-/spin-aware seed &#183; lowest-conformer selection",
         FILL_ML, EDGE_ML, ML),
        ("dft", 3, "DFT geometry + SMD(MeCN), per state [neu/ox1/ox2]",
         "r&#178;SCAN-D4 / def2-SVP(D) &#183; geomeTRIC<BR/>"
         "xTB-ALPB warm start &#183; 3 states concurrent", FILL_GPU, EDGE_GPU, GPU),
        ("sp", 4, "Refined single-point energies",
         "&#969;B97M-V (VV10) / def2-TZVPD &#183; SMD(MeCN)<BR/>"
         "one diffuse basis for all states + Fc", FILL_GPU, EDGE_GPU, GPU),
        ("thermo", 5, "Thermochemistry &#8212; two tiers",
         "pool-wide xTB-RRHO (ranking)<BR/>"
         "native DFT qRRHO (lead confirmation, certificate-locked)",
         FILL_GPU, EDGE_GPU, GPU + " + cert. CPU"),
        ("props", 6, "Derived properties (kept separate)",
         "E&#176; vs Fc &#183; usable capacity &#183; &#955;&#8338; &#183; "
         "&#955;&#7522; (SI)<BR/>disproportionation &#8212; never combined",
         FILL_CPU, EDGE_CPU, CPU),
        ("qc", 7, "QC coverage + cross-point consistency",
         "charge/spin &#183; &#10216;S&#178;&#10217; &#183; geometry &#183; "
         "full state coverage<BR/>failures HELD, not dropped",
         FILL_GATE, EDGE_GATE, CPU),
        ("score", 8, "Scoring / Pareto",
         "&#963;-aware fronts &#183; anolyte &amp; catholyte pools separate<BR/>"
         "<B>no scalar / weighted-sum score</B>", FILL_GATE, EDGE_GATE, CPU),
        ("audit", 9, "Independent table audit",
         "arithmetic re-derivation of exported properties<BR/>"
         "zero-discrepancy acceptance", FILL_GATE, EDGE_GATE, CPU),
    ]
    for name, n, title, detail, fill, edge, tag in stages:
        _stage(g, name, n, title, detail, fill, edge, tag)

    g.edge("title", "build", style="invis")
    order = [s[0] for s in stages]
    edge_labels = {
        ("build", "uma"): " conformer + (q, M)",
        ("uma", "dft"): " gas-phase seed",
        ("dft", "sp"): " solvated geometry",
        ("sp", "thermo"): " E in solvent",
        ("thermo", "props"): " G = E + G_thermal",
        ("props", "qc"): " descriptors",
        ("qc", "score"): " QC-clean couples",
        ("score", "audit"): " exported tables",
    }
    for a, b in zip(order, order[1:]):
        g.edge(a, b, label=edge_labels.get((a, b), ""))

    # ---- inset: level of theory / standard state (manuscript Table 1 digest)
    rows = [
        ("Geometry", "r&#178;SCAN-D4 / def2-SVP(D), SMD(MeCN)"),
        ("Energy", "&#969;B97M-V(VV10) / def2-TZVPD, SMD(MeCN)"),
        ("SCF", "RI-J &#183; conv_tol 1&#215;10&#8315;&#8313; &#183; multi-guess open shell"),
        ("Solvent", "SMD MeCN &#183; &#949;&#8347;=37.5, &#949;&#8338;=1.806, n=1.344"),
        ("Standard state", "298.15 K &#183; RT/F=0.025693 V"),
        ("Reference", "Fc/Fc&#8314;, computed live, level-matched"),
        ("Electrolyte", "TBAPF&#8326; / MeCN &#183; window [&#8722;3.0, 2.0] V"),
        ("Thermo (rank)", "xTB-RRHO on DFT+SMD geometry"),
        ("Thermo (lead)", "native DFT Hessian + projected qRRHO"),
    ]
    trows = "".join(
        f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="18"><B>{k}</B></FONT></TD>'
        f'<TD ALIGN="LEFT"><FONT POINT-SIZE="18" COLOR="{MUTE}">{v}</FONT></TD></TR>'
        for k, v in rows)
    table_label = (
        f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="5" CELLSPACING="0">'
        f'<TR><TD COLSPAN="2"><FONT POINT-SIZE="20"><B>Level of theory &amp; '
        f'conventions</B></FONT></TD></TR>'
        f'{trows}</TABLE>>'
    )
    g.node("methods", label=table_label, shape="box",
           style="rounded,filled,setlinewidth(1.5)", fillcolor="#FFF9E6",
           color="#C89A2B", penwidth="1.8", margin="0.16,0.12")

    # legend
    legend = (
        '<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="6" CELLSPACING="2">'
        '<TR><TD BGCOLOR="' + FILL_GPU + '"><FONT POINT-SIZE="18">&#160;GPU / DFT&#160;</FONT></TD></TR>'
        '<TR><TD BGCOLOR="' + FILL_ML + '"><FONT POINT-SIZE="18">&#160;GPU / ML (UMA)&#160;</FONT></TD></TR>'
        '<TR><TD BGCOLOR="' + FILL_CPU + '"><FONT POINT-SIZE="18">&#160;CPU&#160;</FONT></TD></TR>'
        '<TR><TD BGCOLOR="' + FILL_GATE + '"><FONT POINT-SIZE="18">&#160;QC / scoring&#160;</FONT></TD></TR>'
        '</TABLE>>'
    )
    g.node("legend", label=legend, shape="box",
           style="rounded,filled,setlinewidth(1.2)", fillcolor="white",
           color="#AAAAAA", penwidth="1.2", margin="0.12,0.10")

    # pull the inset to the right of the mid-flow, legend above it
    g.edge("sp", "methods", style="invis", constraint="false")
    g.edge("methods", "legend", style="invis")
    with g.subgraph() as s:
        s.attr(rank="same")
        s.node("dft"); s.node("methods")
    with g.subgraph() as s:
        s.attr(rank="same")
        s.node("uma"); s.node("legend")
    return g


def main():
    if shutil.which("dot") is None:
        print("ERROR: the `dot` executable (Graphviz) is not on PATH.",
              file=sys.stderr)
        sys.exit(1)
    g = build()
    g.attr(dpi="300")
    g.format = "png"  # PNG only, per the manuscript style contract
    g.render(filename=str(FIGDIR / "fig2_pipeline"), cleanup=True)
    print("  ->", FIGDIR / "fig2_pipeline.png")


if __name__ == "__main__":
    main()
