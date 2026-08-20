#!/usr/bin/env python
"""Shared publication style for all pipeline figures.

Import `apply_style()` at the top of every plotting script so the whole figure
set is visually consistent: one font stack, one palette, hairline spines, 300 dpi.
Physics, not fitting — these are cosmetics only; no script here alters a number.
"""
from __future__ import annotations
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DPI = 300

# Palette — colour-blind-safe, muted for print. Keyed by meaning, not by hue.
C = dict(
    dft="#2166AC",       # DFT+SMD (the trusted level)          — blue
    uma="#D6604D",       # UMA gas-phase fallback               — warm red
    measured="#4D4D4D",  # experimental measurement             — dark grey
    good="#1B7837",      # within target / p-type accent        — green
    warn="#B2182B",      # out of band / caution                — red
    accent="#762A83",    # highlight                            — purple
    grid="#B8B8B8",
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
        "savefig.pad_inches": 0.05,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "axes.linewidth": 0.9,
        "axes.edgecolor": "#333333",
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#CCCCCC",
        "mathtext.default": "regular",
    })


def grid_y(ax, alpha=0.35):
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color=C["grid"], lw=0.6, alpha=alpha)


def grid_xy(ax, alpha=0.3):
    ax.set_axisbelow(True)
    ax.grid(True, color=C["grid"], lw=0.6, alpha=alpha)
