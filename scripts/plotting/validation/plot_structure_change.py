#!/usr/bin/env python
"""Structure-change figure (18 pt / 300 dpi, no overlaps; see plot_style.py):

  structure_change.png  heavy-atom RMSD between the two geometries of each electron-
                        transfer event (Kabsch, DFT+SMD geometries). Large bars flag
                        events with big inner-sphere relaxation.

Reads results/structure_descriptors.csv (redox.properties.descriptors) for the candidate
monomers (family != "validation"). The candidate E deg landscape lives in
scripts/plotting/candidates/plot_candidate_analysis.py.

  python scripts/plotting/validation/plot_structure_change.py
"""
from __future__ import annotations
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, C, FAMILY_COLOR, grid_x  # noqa: E402

import matplotlib.pyplot as plt      # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np                   # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = ROOT / "results" / "figures"

# Compact, human-readable names for the candidate monomers (keyed by id).
_PRETTY = {
    "pyridinium":    "N-benzylpyridinium",
    "cyanopyridinium": "N-benzyl-4-cyanopyridinium",
    "viologen":      "Benzyl-methyl viologen",
    "phenothiazine": "N-benzylphenothiazine",
    "tempo":         "4-(benzyloxy)-TEMPO",
    "anthraquinone": "2-(benzyloxymethyl)anthraquinone",
}
# Compact charge-state labels per event (which couple the E deg refers to).
_EVENT = {
    "ox->neu": "+ / 0", "ox2->ox1": "2+ / +", "ox1->neu": "+ / 0",
    "ox->rad": "+ / 0", "rad->red": "0 / -", "neu->red1": "0 / -",
    "red1->red2": "- / 2-", "ox->red": "+ / 0",
}
# Human labels for families (legend order = oxidizing -> reducing character).
_FAMILY_ORDER = [
    "amine (p-type)", "nitroxide", "pyridine", "pyridine-multi-e", "quinone (n-type)",
]


def _read_csv(path: Path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def _fnum(row, key):
    v = (row.get(key) or "").strip()
    if v == "":
        return None
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except ValueError:
        return None


def _candidates(rows):
    return [r for r in rows if r.get("family", "") != "validation"]


def _label(rid, event):
    name = _PRETTY.get(rid, rid)
    couple = _EVENT.get(event, event)
    return f"{name}\n{couple}"


# ---------------------------------------------------------------- structure change
def plot_structure_change(desc_rows, redox_rows):
    fam_by_id = {r["id"]: r.get("family", "") for r in redox_rows}
    pts = []  # (label, rmsd, family)
    for r in desc_rows:
        rid = r["id"]
        if fam_by_id.get(rid, "") == "validation":
            continue
        v = _fnum(r, "rmsd_heavy")
        if v is None:
            continue
        pts.append((_label(rid, r["event"]), v, fam_by_id.get(rid, "")))

    apply_style()
    if not pts:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "No RMSD descriptors available.",
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return fig

    pts.sort(key=lambda t: t[1])  # smallest change at bottom
    labels = [t[0] for t in pts]
    vals = np.array([t[1] for t in pts])
    fams = [t[2] for t in pts]
    colors = [FAMILY_COLOR.get(f, "#888888") for f in fams]

    y = np.arange(len(pts))
    fig, ax = plt.subplots(figsize=(12.5, 0.95 * len(pts) + 2.6))
    ax.barh(y, vals, 0.62, color=colors, edgecolor="k", lw=0.7, zorder=3)

    pad = 0.01 * vals.max()
    for yi, v in zip(y, vals):
        ax.text(v + pad, yi, f"{v:.2f}", ha="left", va="center",
                fontsize=16, color="#111111", zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(pts) - 0.3)
    ax.set_xlabel(r"heavy-atom RMSD between redox states  ($\mathrm{\AA}$)")
    ax.set_title("Structural change on electron transfer (DFT+SMD geometries)")
    grid_x(ax)
    ax.set_xlim(0, vals.max() * 1.16)

    present = [f for f in _FAMILY_ORDER if f in set(fams)]
    handles = [Patch(fc=FAMILY_COLOR[f], ec="k", lw=0.6, label=f) for f in present]
    ax.legend(handles=handles, loc="lower right", title="redox family",
              framealpha=0.95)

    ax.text(0.62, 0.58,
            "large RMSD → big inner-sphere\nrelaxation (single-conformer\nassumption weakest here)",
            transform=ax.transAxes, va="center", ha="left", fontsize=14,
            color="#555555", style="italic")
    return fig


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    redox_rows = _read_csv(ROOT / "results" / "redox_potentials.csv")
    desc_rows = _read_csv(ROOT / "results" / "structure_descriptors.csv")
    n_cand = len(_candidates(redox_rows))
    print(f"redox rows: {len(redox_rows)}  (candidates: {n_cand})   "
          f"descriptor rows: {len(desc_rows)}")

    (FIGDIR / "validation").mkdir(parents=True, exist_ok=True)
    fig = plot_structure_change(desc_rows, redox_rows)
    fig.savefig(FIGDIR / "validation" / "structure_change.png")
    plt.close(fig)
    print(f"  -> {FIGDIR/'validation'/'structure_change.png'}")


if __name__ == "__main__":
    main()
