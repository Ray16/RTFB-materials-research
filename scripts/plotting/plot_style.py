#!/usr/bin/env python
"""Shared publication style for all pipeline figures.

House rules (from the project owner):
  * one panel per figure,
  * 18 pt text everywhere, nothing smaller,
  * no overlapping elements,
  * 300 dpi.

Import `apply_style()` at the top of every plotting script so the whole figure
set is visually consistent. Physics, not fitting — these are cosmetics only;
no script here alters a number.
"""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DPI = 300
BASE = 18  # pt — the floor for every piece of text in the figure set

# Palette — colour-blind-safe, muted for print. Keyed by meaning, not by hue.
C = dict(
    dft="#2166AC",       # DFT+SMD (the trusted level)          — blue
    uma="#D6604D",       # UMA gas-phase fallback               — warm red
    measured="#4D4D4D",  # experimental measurement             — dark grey
    good="#1B7837",      # within target / p-type accent        — green
    warn="#B2182B",      # out of band / caution                — red
    accent="#762A83",    # highlight                            — purple
    grid="#C2C2C2",
    band="#1B7837",
)

# Per-redox-family colours for the candidate landscape (stable, print-friendly).
FAMILY_COLOR = {
    "pyridine":          "#4393C3",
    "pyridine-multi-e":  "#2166AC",
    "amine (p-type)":    "#762A83",
    "nitroxide":         "#E08214",
    "quinone (n-type)":  "#1B7837",
    "validation":        "#4D4D4D",
}

# Okabe-Ito colorblind-safe palette (named), and the redox-family colors used by the
# candidate / D3TaLES reorg figures. Consolidated here so every plot script shares one source.
OKABE_ITO = {
    "blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00",
    "skyblue": "#56B4E9", "reddishpurple": "#CC79A7", "yellow": "#F0E442", "black": "#000000",
}
FAM_COLOR = {"quinone": "#0072B2", "imide": "#009E73", "viologen": "#CC79A7"}


def apply_style():
    plt.rcParams.update({
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.12,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": BASE,
        "axes.titlesize": BASE + 2,
        "axes.titleweight": "bold",
        "axes.labelsize": BASE,
        "axes.labelweight": "medium",
        "axes.linewidth": 1.1,
        "axes.edgecolor": "#333333",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 1.1,
        "ytick.major.width": 1.1,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.labelsize": BASE,
        "ytick.labelsize": BASE,
        "legend.fontsize": BASE,
        "legend.frameon": True,
        "legend.framealpha": 0.95,
        "legend.edgecolor": "#CCCCCC",
        "legend.borderpad": 0.6,
        "legend.labelspacing": 0.5,
        "lines.linewidth": 2.0,
        "patch.linewidth": 0.8,
        "mathtext.default": "regular",
    })


def parity_panel(ax, xs, ys, colors, lim=(0.0, 1.4), band=0.2, point_size=34,
                 marginals=True, marg_color="#9AA3AD"):
    """A polished parity (y-vs-x) panel shared by the reorg figures.

    Fixes overplotting (translucent, edgeless points so the dense core reads as a gradient),
    draws the y=x identity plus a shaded +/-`band` eV agreement envelope, and (optionally) adds
    slim marginal histograms on top and right. Returns a stats dict (n, mad, median, r, frac_in_band).
    Cosmetics only — no number is altered.
    """
    import numpy as np
    lo, hi = lim
    xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    # agreement envelope + identity
    ax.fill_between([lo, hi], [lo - band, hi - band], [lo + band, hi + band],
                    color="#B8C0C8", alpha=0.28, lw=0, zorder=1)
    ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=2)
    ax.scatter(xs, ys, c=colors, s=point_size, alpha=0.62, edgecolors="none", zorder=3)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_axisbelow(True)
    ax.grid(True, color=C["grid"], lw=0.6, alpha=0.30)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # stats
    d = ys - xs
    out = dict(n=len(xs), mad=float(np.mean(np.abs(d))), median=float(np.median(np.abs(d))),
               bias=float(np.mean(d)),
               r=float(np.corrcoef(xs, ys)[0, 1]) if len(xs) > 2 else float("nan"),
               frac_in_band=float(np.mean(np.abs(d) <= band)))
    if marginals:
        bins = np.linspace(lo, hi, 26)
        axt = ax.inset_axes([0, 1.008, 1, 0.15]); axr = ax.inset_axes([1.008, 0, 0.15, 1])
        axt.hist(xs, bins=bins, color=marg_color, alpha=0.85, lw=0)
        axr.hist(ys, bins=bins, orientation="horizontal", color=marg_color, alpha=0.85, lw=0)
        axt.set_xlim(lo, hi); axr.set_ylim(lo, hi)
        for a in (axt, axr):
            a.axis("off")
    return out


def grid_y(ax, alpha=0.4):
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=C["grid"], lw=0.7, alpha=alpha)


def grid_x(ax, alpha=0.4):
    ax.set_axisbelow(True)
    ax.grid(True, axis="x", color=C["grid"], lw=0.7, alpha=alpha)


def grid_xy(ax, alpha=0.35):
    ax.set_axisbelow(True)
    ax.grid(True, color=C["grid"], lw=0.7, alpha=alpha)
