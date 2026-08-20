#!/usr/bin/env python
"""Discussion figures from the redox pipeline outputs.

Reads (all optional-tolerant):
  results/redox_potentials.csv       computed E° table. DFT+SMD columns
                                     (dG_smd_eV, E_vs_Fc_V) may be empty on an early
                                     pass -> we fall back to the UMA gas-phase columns
                                     (dG_gas_uma_eV, E_vs_Fc_gas_V) and label the source.
  results/structure_descriptors.csv  heavy-atom RMSD between adjacent redox states;
                                     inner-sphere lambda if a column is present.
  config/validation.py               experimental anchors (V vs Fc/Fc+, MeCN).

Writes (results/figures/):
  validation_parity.png    computed vs measured E_vs_Fc, y=x guide, MAE/RMSE
  reaction_free_energy.png reaction free energy dG per redox event
  rmsd_by_event.png        heavy-atom RMSD per event (structural change on ET)
  reorg_lambda.png         inner-sphere reorganization energy   (only if data present)
  summary.png              all available panels on one canvas

Physics, not fitting: the parity plot exists to EXPOSE discrepancy with measurement,
not to hide it. We never rescale computed values to match experiment.

  python scripts/plot_results.py
"""
from __future__ import annotations
import csv
import importlib.util
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "results" / "figures"

# --- source-preference: DFT+SMD is the trusted level; UMA gas is a provisional fallback.
E_COLS = [("E_vs_Fc_V", "DFT+SMD"), ("E_vs_Fc_gas_V", "UMA gas")]
DG_COLS = [("dG_smd_eV", "DFT+SMD"), ("dG_gas_uma_eV", "UMA gas")]
# possible names for an inner-sphere reorganization-energy column (schema not yet fixed)
LAMBDA_COLS = ["lambda_inner_eV", "lambda_eV", "lambda_inner", "lambda"]

SRC_COLOR = {"DFT+SMD": "#1f77b4", "UMA gas": "#ff7f0e", "mixed": "#7f7f7f"}


def _load_validation():
    spec = importlib.util.spec_from_file_location(
        "validation", ROOT / "config" / "validation.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.VALIDATION


def _read_csv(path: Path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def _fnum(row, key):
    """Float value or None for an empty/missing/NaN cell."""
    v = (row.get(key) or "").strip()
    if v == "":
        return None
    try:
        x = float(v)
        return None if math.isnan(x) else x
    except ValueError:
        return None


def _pick(row, cols):
    """First non-empty value across a preference list of (col, source_label)."""
    for col, src in cols:
        v = _fnum(row, col)
        if v is not None:
            return v, src
    return None, None


def _label(rid, event):
    return f"{rid}\n{event}"


# ---------------------------------------------------------------- parity plot
def plot_parity(redox_rows, validation, ax=None):
    """Computed vs measured E_vs_Fc. Match by (id, event). y=x guide + MAE/RMSE.

    Only validation-set ids that also appear in the computed table are plotted; the
    decorated monomers have no experimental anchor and are intentionally excluded here.
    """
    by_key = {(r["id"], r["event"]): r for r in redox_rows}
    pts = []  # (exp, comp, src, label)
    for entry in validation:
        rid = entry["id"]
        for ev in entry.get("events", []):
            key = (rid, ev["event"])
            row = by_key.get(key)
            if row is None:
                continue
            comp, src = _pick(row, E_COLS)
            if comp is None:
                continue
            pts.append((ev["exp_V_vs_Fc"], comp, src, _label(rid, ev["event"])))

    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(6.4, 6.0))

    if not pts:
        ax.text(0.5, 0.5, "No validation anchors matched yet.\n"
                          "(Re-run redox.py over the validation cores\n"
                          "once their DFT+SMD finishes.)",
                ha="center", va="center", transform=ax.transAxes, fontsize=10)
        ax.set_title("Validation parity — awaiting data")
        return (fig if own else None), pts

    exp = np.array([p[0] for p in pts])
    comp = np.array([p[1] for p in pts])
    srcs = [p[2] for p in pts]

    lo = min(exp.min(), comp.min()) - 0.3
    hi = max(exp.max(), comp.max()) + 0.3
    ax.plot([lo, hi], [lo, hi], "--", color="0.5", lw=1, zorder=1, label="y = x")

    for src in sorted(set(srcs)):
        m = [i for i, s in enumerate(srcs) if s == src]
        ax.scatter(exp[m], comp[m], s=70, color=SRC_COLOR.get(src, "0.3"),
                   edgecolor="k", lw=0.5, zorder=3, label=src)
    for e, c, s, lab in pts:
        ax.annotate(lab.replace("\n", " "), (e, c), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")

    err = comp - exp
    mae = float(np.mean(np.abs(err)))
    rmse = float(np.sqrt(np.mean(err ** 2)))
    src_note = "+".join(sorted(set(srcs)))
    ax.text(0.03, 0.97, f"n = {len(pts)}   source: {src_note}\n"
                        f"MAE  = {mae:.3f} V\nRMSE = {rmse:.3f} V",
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Experimental  E  (V vs Fc/Fc$^+$)")
    ax.set_ylabel("Computed  E  (V vs Fc/Fc$^+$)")
    ax.set_title("Validation parity: computed vs measured")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    return (fig if own else None), pts


# ---------------------------------------------------------- bar: reaction dG
def plot_dg(redox_rows, ax=None):
    rows = []
    for r in redox_rows:
        v, src = _pick(r, DG_COLS)
        if v is not None:
            rows.append((_label(r["id"], r["event"]), v, src))
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(rows) + 2), 4.6))
    if not rows:
        ax.text(0.5, 0.5, "No reaction free energies available.",
                ha="center", va="center", transform=ax.transAxes)
        return fig if own else None
    labels = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    cols = [SRC_COLOR.get(r[2], "0.3") for r in rows]
    x = np.arange(len(rows))
    ax.bar(x, vals, color=cols, edgecolor="k", lw=0.4)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(r"$\Delta G$ of reduction  (eV)")
    srcset = sorted({r[2] for r in rows})
    ax.set_title(f"Reaction free energy per redox event  (source: {'+'.join(srcset)})")
    ax.grid(True, axis="y", alpha=0.3)
    handles = [plt.Rectangle((0, 0), 1, 1, color=SRC_COLOR[s]) for s in srcset]
    ax.legend(handles, srcset, fontsize=8)
    return fig if own else None


# ------------------------------------------------------------- bar: RMSD
def plot_rmsd(desc_rows, ax=None):
    rows = [(_label(r["id"], r["event"]), _fnum(r, "rmsd_heavy"),
             r.get("geom_source", "")) for r in desc_rows
            if _fnum(r, "rmsd_heavy") is not None]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(rows) + 2), 4.6))
    if not rows:
        ax.text(0.5, 0.5, "No RMSD descriptors available.",
                ha="center", va="center", transform=ax.transAxes)
        return fig if own else None
    gsrc_color = {"dft_smd": "#1f77b4", "uma_gas": "#ff7f0e"}
    labels = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    cols = [gsrc_color.get(r[2], "0.5") for r in rows]
    x = np.arange(len(rows))
    ax.bar(x, vals, color=cols, edgecolor="k", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(r"heavy-atom RMSD  ($\AA$)")
    gset = sorted({r[2] for r in rows if r[2]})
    ax.set_title("Structural change on electron transfer (Kabsch heavy-atom RMSD)")
    ax.grid(True, axis="y", alpha=0.3)
    if gset:
        handles = [plt.Rectangle((0, 0), 1, 1, color=gsrc_color.get(s, "0.5"))
                   for s in gset]
        ax.legend(handles, gset, fontsize=8, title="geometry")
    return fig if own else None


# ------------------------------------------------------- bar: lambda (optional)
def _lambda_col(desc_rows):
    if not desc_rows:
        return None
    cols = desc_rows[0].keys()
    for c in LAMBDA_COLS:
        if c in cols:
            return c
    return None


def plot_lambda(desc_rows, ax=None):
    col = _lambda_col(desc_rows)
    rows = []
    if col:
        rows = [(_label(r["id"], r["event"]), _fnum(r, col)) for r in desc_rows
                if _fnum(r, col) is not None]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(max(7, 0.7 * max(len(rows), 1) + 2), 4.6))
    if not rows:
        ax.text(0.5, 0.5, "Inner-sphere reorganization energy not yet computed.\n"
                          "(descriptors.lambda_inner is implemented but not wired\n"
                          "into the CSV — enable it once DFT geometries exist.)",
                ha="center", va="center", transform=ax.transAxes, fontsize=9)
        ax.set_title(r"Reorganization energy $\lambda$ — awaiting data")
        return fig if own else None
    labels = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    x = np.arange(len(rows))
    ax.bar(x, vals, color="#2ca02c", edgecolor="k", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel(r"$\lambda_{inner}$  (eV)")
    ax.set_title(r"Inner-sphere reorganization energy $\lambda$")
    ax.grid(True, axis="y", alpha=0.3)
    return fig if own else None


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    redox_rows = _read_csv(ROOT / "results" / "redox_potentials.csv")
    desc_rows = _read_csv(ROOT / "results" / "structure_descriptors.csv")
    validation = _load_validation()

    print(f"redox rows: {len(redox_rows)}   descriptor rows: {len(desc_rows)}")

    # individual figures
    fig, pts = plot_parity(redox_rows, validation)
    fig.tight_layout(); fig.savefig(FIGDIR / "validation_parity.png", dpi=150)
    plt.close(fig)
    print(f"validation_parity.png  ({len(pts)} anchors matched)")

    for name, fn, src in [
        ("reaction_free_energy.png", plot_dg, redox_rows),
        ("rmsd_by_event.png", plot_rmsd, desc_rows),
        ("reorg_lambda.png", plot_lambda, desc_rows),
    ]:
        f = fn(src)
        f.tight_layout(); f.savefig(FIGDIR / name, dpi=150); plt.close(f)
        print(name)

    # combined summary canvas
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(2, 2)
    plot_parity(redox_rows, validation, ax=fig.add_subplot(gs[0, 0]))
    plot_dg(redox_rows, ax=fig.add_subplot(gs[0, 1]))
    plot_rmsd(desc_rows, ax=fig.add_subplot(gs[1, 0]))
    plot_lambda(desc_rows, ax=fig.add_subplot(gs[1, 1]))
    fig.suptitle("Redox screening — pipeline results overview", fontsize=14, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(FIGDIR / "summary.png", dpi=150); plt.close(fig)
    print("summary.png")
    print(f"\nFigures written to {FIGDIR}")


if __name__ == "__main__":
    main()
