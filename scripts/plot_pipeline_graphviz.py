#!/usr/bin/env python
"""Workflow diagram (Graphviz edition): how a molecule flows through the pipeline.

Rebuilds the old hand-laid-matplotlib flowchart (scripts/plot_pipeline.py) with
Graphviz, which does the layout for us — cleaner routing, no hand-tuned
coordinates, and crisp records. The two inputs (candidate monomers + the known-E
validation set) share the IDENTICAL pipeline; the validation gate is what
licenses the final ranking.

Needs the system `dot` executable (package: graphviz). If it is missing this
script says so and exits non-zero rather than writing a broken file.

  python scripts/plot_pipeline_graphviz.py   ->  results/figures/pipeline.png

Physics, not fitting: this only draws the workflow; it changes no number.
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import C  # noqa: E402  (reuse the house palette)

import graphviz  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "results" / "figures"

# Palette (house style) -> Graphviz fills. Keep dark text on light fills.
INK = "#111111"
MUTE = "#555555"
FILL_INPUT_C = "#EDEDED"
FILL_INPUT_V = "#E6E0EF"
FILL_RDKIT = "#DCE6F1"
FILL_UMA = "#FBE3DE"
FILL_DFT = "#CBD9EC"
FILL_REDOX = "#DCEBDD"
FILL_GATE = "#FDF0D5"
FILL_NOTE = "#FFF9E6"


def _node(g, name, title, detail, fill, edge):
    """A two-line box: bold title over a smaller detail line (HTML label)."""
    label = (
        f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2">'
        f'<TR><TD><FONT POINT-SIZE="15"><B>{title}</B></FONT></TD></TR>'
        f'<TR><TD><FONT POINT-SIZE="11" COLOR="{MUTE}">{detail}</FONT></TD></TR>'
        f'</TABLE>>'
    )
    g.node(name, label=label, fillcolor=fill, color=edge)


def build():
    g = graphviz.Digraph("redox_pipeline", format="png")
    g.attr(rankdir="TB", bgcolor="white", nodesep="0.35", ranksep="0.55",
           splines="true", pad="0.3")
    g.attr("node", shape="box", style="rounded,filled,setlinewidth(2)",
           fontname="DejaVu Sans", fontcolor=INK, penwidth="2",
           margin="0.18,0.10")
    g.attr("edge", color="#3A3A3A", penwidth="2", arrowsize="0.9",
           fontname="DejaVu Sans", fontsize="11", fontcolor="#333333")

    # title (a plain node up top)
    g.node("title",
           label=('<<FONT POINT-SIZE="20"><B>Redox-potential screening '
                  'pipeline</B></FONT><BR/>'
                  '<FONT POINT-SIZE="12" COLOR="' + MUTE + '">'
                  'ML gas-phase pre-opt &#8594; DFT geometry &amp; energy in '
                  'implicit solvent &#8594; E&#176; vs Fc/Fc&#8314;</FONT>>'),
           shape="plaintext", style="", margin="0")

    # ---- two inputs, same rank, sharing one pipeline
    _node(g, "cand", "Candidate monomers", "benzylic-site redox groups",
          FILL_INPUT_C, "#7A7A7A")
    _node(g, "val", "Validation set", "cores with measured E&#176;",
          FILL_INPUT_V, "#8C6BA8")
    with g.subgraph() as s:
        s.attr(rank="same")
        s.node("cand"); s.node("val")

    # charge/spin note (diamond-free; a highlighted box on the merge)
    g.node("cs",
           label=('<<FONT POINT-SIZE="12"><I>each redox state = '
                  '(charge <B>q</B>, spin multiplicity)<BR/>'
                  '&#8212; never assume neutral singlet</I></FONT>>'),
           shape="note", style="filled", fillcolor=FILL_NOTE,
           color="#E0C36A", fontcolor="#333333", penwidth="1.5")

    # ---- main flow
    _node(g, "rdkit", "3D structure generation",
          "RDKit ETKDG conformer ensemble &#183; stereo verified",
          FILL_RDKIT, "#4B7BB0")
    _node(g, "uma", "UMA ML pre-optimization (gas phase)",
          "uma-s-1p2p1 &#183; per (charge, spin) &#183; GPU-batched",
          FILL_UMA, C["uma"])
    _node(g, "dftopt", "DFT geometry optimization in SMD",
          "r2SCAN-D4 / def2-SVP(D) &#183; MeCN (&#949;=37.5)",
          FILL_DFT, C["dft"])
    _node(g, "dftsp", "DFT single-point energy in SMD",
          "&#969;B97M-V / def2-TZVP(D) + VV10 &#183; diffuse basis for anions",
          FILL_DFT, C["dft"])
    _node(g, "redox", "Redox potential",
          "E&#176; = &#8722;&#916;G/nF &#8722; E(Fc) &#183; vs Fc/Fc&#8314;",
          FILL_REDOX, C["good"])
    _node(g, "gate", "Validation gate: computed vs measured E&#176;",
          "MAE 0.18 V &#183; within &#177;0.15 V band &#8594; trust the ranking",
          FILL_GATE, "#C89A2B")

    # edges
    g.edge("title", "cand", style="invis")
    g.edge("title", "val", style="invis")
    g.edge("cand", "cs")
    g.edge("val", "cs")
    g.edge("cs", "rdkit")
    g.edge("rdkit", "uma", label=" conformer + charge/spin")
    g.edge("uma", "dftopt", label=" gas-phase geometry")
    g.edge("dftopt", "dftsp", label=" solvated geometry")
    g.edge("dftsp", "redox", label=" G in solvent")
    g.edge("redox", "gate")

    # footer
    g.node("foot",
           label=('<<FONT POINT-SIZE="11" COLOR="#666666"><I>No fitting: '
                  'computed potentials are never rescaled to experiment; '
                  'the gate exposes error, it does not hide it.</I></FONT>>'),
           shape="plaintext", style="", margin="0")
    g.edge("gate", "foot", style="invis")
    return g


def main():
    if shutil.which("dot") is None:
        print("ERROR: the `dot` executable (Graphviz) is not on PATH.\n"
              "  Install it (apt: graphviz | conda: graphviz) and re-run.",
              file=sys.stderr)
        sys.exit(1)

    FIGDIR.mkdir(parents=True, exist_ok=True)
    g = build()
    # High DPI for print parity with the rest of the figure set.
    g.attr(dpi="300")
    out = FIGDIR / "pipeline"          # graphviz appends .png
    g.render(filename=str(out), cleanup=True)
    print(f"  -> {out.with_suffix('.png')}")


if __name__ == "__main__":
    main()
