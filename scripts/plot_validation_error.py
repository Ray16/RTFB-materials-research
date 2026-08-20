#!/usr/bin/env python
"""Validation figure (hard gate): measured vs computed E° per experimental redox event.

Single panel, horizontal grouped bars — one pair per validation event: the measured
value (dark grey) beside our computed value (family/level colour). Because the bars are
horizontal, the event labels sit on the y-axis and never rotate or overlap. The signed
residual is printed at the end of each pair and an MAE/RMSE box states the gate outcome.

We NEVER rescale to match experiment — the bars show the true residual, warts and all.

  results/redox_potentials.csv   computed E° table (needs E_vs_Fc_V populated)
  config/validation.py           experimental anchors (V vs Fc/Fc+, MeCN)

  python scripts/plot_validation_error.py   ->  results/figures/validation.png
"""
from __future__ import annotations
import csv, importlib.util, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, C, grid_x  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

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


# Pretty, compact labels for the y-axis (avoid raw ids/events colliding).
_PRETTY = {
    "ferrocene":            "Ferrocene",
    "methyl_viologen":      "Methyl viologen",
    "tempo_parent":         "TEMPO",
    "phenothiazine_parent": "Phenothiazine",
    "anthraquinone_parent": "Anthraquinone",
    "methylpyridinium":     "Methylpyridinium",
}
_EVENT = {
    "ox->neu": "ox / neu", "ox2->ox1": "2+ / +", "ox1->neu": "+ / 0",
    "ox->rad": "ox / rad", "neu->red1": "0 / -", "red1->red2": "- / 2-",
    "ox->red": "ox / red",
}


def main():
    apply_style()
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
            name = _PRETTY.get(rid, rid)
            evl = _EVENT.get(ev["event"], ev["event"])
            lab = f"{name}\n{evl}"
            pts.append((lab, exp, comp, err, src))
            print(f"  {rid:22s} {ev['event']:12s} {exp:+7.3f} {comp:+7.3f} {err:+7.3f}  {src}")

    if not pts:
        fig, ax = plt.subplots(figsize=(9, 6))
        ax.text(0.5, 0.5, "No validation events scored yet.",
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        fig.savefig(FIGDIR / "validation.png")
        print("\n  no anchors scored — placeholder written"); return

    # order by measured potential (most reducing at bottom -> most oxidizing at top)
    pts.sort(key=lambda t: t[1])
    labels = [t[0] for t in pts]
    exps = np.array([t[1] for t in pts])
    comps = np.array([t[2] for t in pts])
    errs = np.array([t[3] for t in pts])
    srcs = [t[4] for t in pts]

    mae = float(np.mean(np.abs(errs))); rmse = float(np.sqrt(np.mean(errs**2)))
    mbe = float(np.mean(errs))
    n_in = int(np.sum(np.abs(errs) <= TARGET_BAND))

    y = np.arange(len(pts))
    h = 0.38
    comp_colors = [C["dft"] if s == "DFT+SMD" else C["uma"] for s in srcs]

    fig, ax = plt.subplots(figsize=(12.5, 1.05 * len(pts) + 3.2))

    ax.barh(y + h / 2, exps, h, color=C["measured"], edgecolor="k", lw=0.6,
            zorder=3, label="Measured (experiment)")
    ax.barh(y - h / 2, comps, h, color=comp_colors, edgecolor="k", lw=0.6,
            zorder=3, label="Computed (DFT+SMD)")

    # value annotations at the tip of each bar
    span = max(abs(comps.min()), abs(comps.max()), abs(exps.min()), abs(exps.max()))
    pad = 0.03 * span
    for yi, v in zip(y + h / 2, exps):
        ax.text(v + (pad if v >= 0 else -pad), yi, f"{v:+.2f}",
                ha="left" if v >= 0 else "right", va="center",
                fontsize=15, color=C["measured"], zorder=4)
    for yi, v, er in zip(y - h / 2, comps, errs):
        ax.text(v + (pad if v >= 0 else -pad), yi, f"{v:+.2f}",
                ha="left" if v >= 0 else "right", va="center",
                fontsize=15, color="#111111", zorder=4)

    ax.axvline(0, color="k", lw=1.1, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.8, len(pts) - 0.2)
    ax.set_xlabel(r"$E^\circ$  (V vs Fc/Fc$^+$)")
    ax.set_title("Validation gate: computed vs measured redox potentials")
    grid_x(ax)

    # widen x so annotations fit
    xlo = min(0, comps.min(), exps.min()) - 0.55
    xhi = max(0, comps.max(), exps.max()) + 0.55
    ax.set_xlim(xlo, xhi)

    ax.legend(loc="lower right", framealpha=0.95)
    ax.text(0.015, 0.985,
            f"n = {len(pts)}   MAE = {mae:.2f} V   RMSE = {rmse:.2f} V\n"
            f"mean signed = {mbe:+.2f} V   {n_in}/{len(pts)} within $\\pm${TARGET_BAND:.2f} V\n"
            "No fitting — residuals shown as computed.",
            transform=ax.transAxes, va="top", ha="left", fontsize=15,
            bbox=dict(boxstyle="round,pad=0.5", fc="#FBFBF7", ec="#B9B9B9", lw=1.2))

    fig.savefig(FIGDIR / "validation.png")
    print(f"\n  MAE = {mae:.3f} V   RMSE = {rmse:.3f} V   mean signed = {mbe:+.3f} V")
    print(f"  -> {FIGDIR/'validation.png'}")


if __name__ == "__main__":
    main()
