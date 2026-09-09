"""Shared publication-ready matplotlib style for all project figures.

Call `set_pub_style()` at the top of every plotting script so figures are consistent and meet
the project standard: >=18 pt text, no overlapping text (constrained layout), clean top-journal
theme (sans-serif, no top/right spines, colorblind-safe Okabe-Ito palette, 300 dpi).
"""
from __future__ import annotations

# Okabe-Ito colorblind-safe palette
OKABE_ITO = {
    "blue": "#0072B2", "orange": "#E69F00", "green": "#009E73", "vermillion": "#D55E00",
    "skyblue": "#56B4E9", "reddishpurple": "#CC79A7", "yellow": "#F0E442", "black": "#000000",
}
# project family colors
FAM_COLOR = {"quinone": "#0072B2", "imide": "#009E73", "viologen": "#CC79A7"}


def set_pub_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        # fonts — everything >= 18 pt
        "font.size": 18,
        "axes.titlesize": 18,
        "axes.labelsize": 18,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 16,
        "figure.titlesize": 20,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "dejavusans",
        # clean top-journal theme
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 1.3,
        "xtick.major.width": 1.3,
        "ytick.major.width": 1.3,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.grid": False,
        "legend.frameon": False,
        # layout — prevent overlaps, export at print quality
        "figure.constrained_layout.use": True,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })
