#!/usr/bin/env python
"""Validation-error bar plot: computed - measured E° per experimental redox event.

This is the hard-gate figure. One signed bar per validation event, referenced to
ferrocene at OUR level (redox.py fills E_vs_Fc_V using electrolyte.FC_ABS_COMPUTED_V).
We NEVER rescale to match experiment — the bars show the true residual, warts and all.

  results/redox_potentials.csv   computed E° table (needs E_vs_Fc_V populated)
  config/validation.py           experimental anchors (V vs Fc/Fc+, MeCN)

Writes results/figures/validation_error_bars.png and prints the numeric table.

  python scripts/plot_validation_error.py
"""
from __future__ import annotations
import csv, importlib.util, math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "results" / "figures"
# DFT+SMD is the trusted level; fall back to UMA gas only if DFT is missing, and label it.
E_COLS = [("E_vs_Fc_V", "DFT+SMD"), ("E_vs_Fc_gas_V", "UMA gas")]
TARGET_BAND = 0.15  # V — typical MAE band for this protocol class; a visual guide, not a fit


def _load_validation():
    s = importlib.util.spec_from_file_location("validation", ROOT / "config" / "validation.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m.VALIDATION


def _fnum(row, key):
    v = (row.get(key) or "").strip()
    if not v:
        return None
    try:
        x = float(v); return None if math.isnan(x) else x
    except ValueError:
        return None


def _pick(row):
    for col, src in E_COLS:
        v = _fnum(row, col)
        if v is not None:
            return v, src
    return None, None


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    p = ROOT / "results" / "redox_potentials.csv"
    rows = list(csv.DictReader(p.open())) if p.exists() else []
    by = {(r["id"], r["event"]): r for r in rows}
    val = _load_validation()

    pts = []  # (label, exp, comp, err, src)
    print(f"\n  {'id':22s} {'event':12s} {'exp':>7s} {'comp':>7s} {'error':>7s}  src")
    for e in val:
        rid = e["id"]
        for ev in e.get("events", []):
            row = by.get((rid, ev["event"]))
            if not row:
                continue
            comp, src = _pick(row)
            if comp is None:
                continue
            exp = ev["exp_V_vs_Fc"]; err = comp - exp
            lab = f"{rid}\n{ev['event']}"
            pts.append((lab, exp, comp, err, src))
            print(f"  {rid:22s} {ev['event']:12s} {exp:+7.3f} {comp:+7.3f} {err:+7.3f}  {src}")

    if not pts:
        fig, ax = plt.subplots(figsize=(7, 5.2))
        ax.text(0.5, 0.5, "No validation events scored yet.\n"
                          "(Run redox.py after the validation-core DFT+SMD finishes.)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_title("Validation error — awaiting data")
        fig.tight_layout(); fig.savefig(FIGDIR / "validation_error_bars.png", dpi=150)
        print("\n  no anchors scored — placeholder written"); return

    labels = [p[0] for p in pts]
    exps = np.array([p[1] for p in pts])
    comps = np.array([p[2] for p in pts])
    errs = np.array([p[3] for p in pts])
    srcs = [p[4] for p in pts]
    x = np.arange(len(pts))

    mae = float(np.mean(np.abs(errs))); rmse = float(np.sqrt(np.mean(errs**2)))
    mbe = float(np.mean(errs))
    srcset = "+".join(sorted(set(srcs)))

    # Two panels: (top) measured vs computed E° side-by-side; (bottom) signed residual.
    fig, (axv, axe) = plt.subplots(
        2, 1, figsize=(max(8, 1.15 * len(pts) + 2), 8.4),
        gridspec_kw=dict(height_ratios=[1.5, 1.0], hspace=0.08), sharex=True)

    # --- top: grouped bars, experimental measurement alongside our computed value ---
    w = 0.4
    axv.bar(x - w / 2, exps, w, color="#555555", edgecolor="k", lw=0.5, zorder=3,
            label="measured (exp, V vs Fc/Fc$^+$)")
    comp_colors = ["#1f77b4" if s == "DFT+SMD" else "#ff7f0e" for s in srcs]
    axv.bar(x + w / 2, comps, w, color=comp_colors, edgecolor="k", lw=0.5, zorder=3,
            label="computed (DFT+SMD)")
    for xi, v in zip(x - w / 2, exps):
        axv.annotate(f"{v:+.2f}", (xi, v), ha="center",
                     va="bottom" if v >= 0 else "top", fontsize=7,
                     xytext=(0, 2 if v >= 0 else -2), textcoords="offset points")
    for xi, v in zip(x + w / 2, comps):
        axv.annotate(f"{v:+.2f}", (xi, v), ha="center",
                     va="bottom" if v >= 0 else "top", fontsize=7,
                     xytext=(0, 2 if v >= 0 else -2), textcoords="offset points")
    axv.axhline(0, color="k", lw=1)
    axv.set_ylabel(r"$E$  (V vs Fc/Fc$^+$)")
    axv.set_title("Validation: measured vs computed redox potential")
    axv.grid(True, axis="y", alpha=0.3)
    axv.legend(loc="best", fontsize=8)
    axv.text(0.02, 0.97,
             f"n = {len(pts)}   source: {srcset}\nMAE = {mae:.3f} V   RMSE = {rmse:.3f} V\n"
             f"mean signed err = {mbe:+.3f} V",
             transform=axv.transAxes, va="top", ha="left", fontsize=9,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    # --- bottom: signed residual (computed - measured) ---
    axe.bar(x, errs, color=comp_colors, edgecolor="k", lw=0.5, zorder=3)
    axe.axhline(0, color="k", lw=1)
    axe.axhspan(-TARGET_BAND, TARGET_BAND, color="green", alpha=0.10, zorder=0,
                label=f"±{TARGET_BAND:.2f} V target band")
    for xi, ei in zip(x, errs):
        axe.annotate(f"{ei:+.2f}", (xi, ei), ha="center",
                     va="bottom" if ei >= 0 else "top", fontsize=8,
                     xytext=(0, 3 if ei >= 0 else -3), textcoords="offset points")
    axe.set_xticks(x); axe.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    axe.set_ylabel(r"computed $-$ measured (V)")
    axe.set_title("Residual per redox event  (positive = computed too oxidizing)")
    axe.grid(True, axis="y", alpha=0.3)
    axe.legend(loc="lower right", fontsize=8)

    fig.tight_layout(); fig.savefig(FIGDIR / "validation_error_bars.png", dpi=150)
    print(f"\n  MAE = {mae:.3f} V   RMSE = {rmse:.3f} V   mean signed = {mbe:+.3f} V")
    print(f"  -> {FIGDIR/'validation_error_bars.png'}")


if __name__ == "__main__":
    main()
