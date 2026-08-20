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


def grid_y(ax, alpha=0.4):
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=C["grid"], lw=0.7, alpha=alpha)


def grid_x(ax, alpha=0.4):
    ax.set_axisbelow(True)
    ax.grid(True, axis="x", color=C["grid"], lw=0.7, alpha=alpha)


def grid_xy(ax, alpha=0.35):
    ax.set_axisbelow(True)
    ax.grid(True, color=C["grid"], lw=0.7, alpha=alpha)
