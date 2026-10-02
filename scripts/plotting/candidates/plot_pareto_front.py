#!/usr/bin/env python
"""Visualize the sigma-aware multi-objective Pareto front (redox.screening.pareto), per capacity/salt
scenario (config CAPACITY_SCENARIOS). Three complementary views, PNG only:

  1. pareto_parallel_coords.png   one vertical axis per objective (up = better), one line per
                                  candidate; axes min-max scaled FOR DISPLAY ONLY (raw range
                                  printed at each axis end); sigma whiskers on front members.
  2. pareto_scatter_matrix_<scenario>.png   every objective pair with sigma error bars.
  3. pareto_dominance_heatmap.png candidate x objective, coloured by within-pool rank, raw
                                  values annotated, plus who dominates each off-front candidate.

Encoding shared by all three: each candidate has ONE fixed colour (Okabe-Ito, colour-blind
safe); front = filled marker / solid line, dominated = hollow marker / dashed line. A single
legend replaces per-point text labels (no overlaps).

Objectives (higher = better after orientation): voltage = -E (anolyte), capacity (scenario
column), stability = dG_disp, kinetics = -lambda_i. lambda_o / lambda_het / lambda_se are annotations,
not objectives. No scalarized figure of merit.

  python scripts/plotting/candidates/plot_pareto_front.py
"""
from __future__ import annotations
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
import sys  # noqa: E402
sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
sys.path.insert(0, str(ROOT / "src"))
import plot_style as ps  # noqa: E402
from redox.core.common import load_config as _cfg  # noqa: E402

OUT = ROOT / "results" / "figures" / "candidates"
POOL = "anolyte"
# (key, label, unit, sign): display value = sign * raw, so up/right = better
AXES = [("voltage", "voltage  −E", "V vs Fc", -1.0),
        ("capacity", "capacity", "mAh g$^{-1}$", 1.0),
        ("stability", "stability  ΔG$_{disp}$", "kJ mol$^{-1}$", 1.0),
        ("kinetics", "kinetics  −λ$_i$", "eV", -1.0)]
SIGN = {a[0]: a[3] for a in AXES}
OKABE_ITO = ["#D55E00", "#009E73", "#0072B2", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
FS_TICK, FS_LABEL, FS_TITLE, FS_LEG, FS_CELL = 16, 18, 20, 16, 16


def _f(x):
    try:
        return float(x) if x not in (None, "") else None
    except ValueError:
        return None


def load():
    sc = _cfg("scorecard_config")
    card = {r["id"]: r for r in csv.DictReader(open(ROOT / "results/scorecard.csv"))}
    data = {}
    for r in csv.DictReader(open(ROOT / "results/pareto_shortlist.csv")):
        if r["pool"] != POOL:
            continue
        c = card[r["id"]]
        raw = dict(voltage=_f(r["E_V"]), capacity=_f(r["capacity"]),
                   stability=_f(r["dG_disp_kJmol"]), kinetics=_f(r["lambda_i_eV"]))
        sig = dict(voltage=_f(c["sigma_E_V"]) or 0.0, capacity=0.0,
                   stability=(_f(c["sigma_disp_eV"]) or sc.SIGMA_DISP_EV) * 96.485,
                   kinetics=_f(c["sigma_lambda_eV"]) or sc.SIGMA_LAMBDA_EV)
        data.setdefault(r["scenario"], []).append(dict(
            id=r["id"], raw=raw, sig=sig, front=(r["pareto_optimal"] == "True"),
            incomplete=bool(r.get("missing_objectives")),
            dominated_by=[n for n in r["dominated_by"].split(";") if n],
            lam_qc=r["lambda_qc"]))
    return data


# One DISTINCT colour per candidate (Okabe-Ito + Paul Tol "muted" extension, both CVD-safe),
# hue loosely by family: quinones = blues/indigo/teal/pink, imides = oranges, viologens =
# green/black/olive. Every current candidate is listed, so none falls through to a shared grey.
FIXED_COLOR = {
    "mophquinone": "#0072B2", "dmophquinone": "#56B4E9", "aq_benzyloxy": "#332288",
    "aq_benzylamino": "#CC79A7", "nq_benzyloxy": "#44AA99",
    "pmdi": "#E69F00", "ndi_ammonium": "#D55E00",
    "viologen": "#009E73", "ethylviologen": "#000000", "bisviologen": "#999933",
    "dtbc_phenol": "#882255",
}
SPARE = ["#AA4499", "#117733", "#DDCC77", "#88CCEE"]


def _palette(data):
    ids = sorted({r["id"] for rows in data.values() for r in rows})
    spare = iter(SPARE)
    out = {i: FIXED_COLOR.get(i) or next(spare, None) for i in ids}
    missing = [i for i, c in out.items() if c is None]
    if missing:
        raise SystemExit(f"plot_pareto_front: no distinct colour left for {missing} — extend FIXED_COLOR")
    return out


def _legend_handles(ids, cmap):
    h = [Line2D([0], [0], color=cmap[i], lw=3, marker="o", ms=11, label=i) for i in ids]
    h += [Line2D([0], [0], color="#333", lw=3, marker="o", ms=11, label="front (filled, solid)"),
          Line2D([0], [0], color="#333", lw=2, ls="--", marker="o", ms=11, mfc="white",
                 label="dominated (hollow, dashed)"),
          Line2D([0], [0], color="#333", lw=2.2, ls=":", marker="s", ms=10, mfc="white",
                 label="incomplete: missing axis (dotted)")]
    return h


def _style_axes(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=FS_TICK)


# ---------------------------------------------------------------- 1. parallel coordinates
def parallel_coords(data, cmap):
    scen = list(data)
    fig, axs = plt.subplots(1, len(scen), figsize=(10.5 * len(scen) + 4, 8.2))
    axs = np.atleast_1d(axs)
    for ax, s in zip(axs, scen):
        rows = data[s]
        x = np.arange(len(AXES))
        lims = {}
        for key, *_ in AXES:
            disp = [SIGN[key] * r["raw"][key] for r in rows if r["raw"][key] is not None]
            lo, hi = min(disp), max(disp)
            pad = 0.06 * (hi - lo if hi > lo else 1.0)
            lims[key] = (lo - pad, hi + pad)

        def norm(key, raw):
            lo, hi = lims[key]
            return (SIGN[key] * raw - lo) / (hi - lo)

        front_ids = sorted(r["id"] for r in rows if r["front"])
        offs = {i: 0.13 * (k - (len(front_ids) - 1) / 2) for k, i in enumerate(front_ids)}
        for r in sorted(rows, key=lambda r: r["front"]):          # front drawn on top
            col = cmap[r["id"]]
            ys = [norm(k, r["raw"][k]) if r["raw"][k] is not None else np.nan
                  for k, *_ in AXES]
            if r["front"]:
                ax.plot(x, ys, "-o", color=col, lw=3.6, ms=12, zorder=4)
                for xi, (k, *_rest) in zip(x, AXES):
                    if r["raw"][k] is None or not r["sig"][k]:
                        continue
                    lo, hi = lims[k]
                    y = norm(k, r["raw"][k]); e = r["sig"][k] / (hi - lo)
                    y0, y1 = max(y - e, -0.02), min(y + e, 1.02)
                    xw = xi + offs[r["id"]]
                    ax.plot([xw, xw], [y0, y1], color=col, lw=2.0, zorder=3)
                    for yy, clipped, mk in ((y0, y - e < -0.02, "v"), (y1, y + e > 1.02, "^")):
                        ax.plot(xw, yy, marker=mk if clipped else "_", color=col,
                                ms=9 if clipped else 12, mew=2, zorder=3)
            elif r["incomplete"]:           # missing a primary axis: NOT ranked, not "dominated"
                ax.plot(x, ys, ":s", color=col, lw=2.2, ms=10, mfc="white", mew=2.0,
                        alpha=0.9, zorder=2)
            else:
                ax.plot(x, ys, "--o", color=col, lw=2.0, ms=10, mfc="white", mew=2.0,
                        alpha=0.85, zorder=2)
        for xi, (k, lab, unit, sign) in zip(x, AXES):
            ax.axvline(xi, color="#444", lw=1.2, zorder=1)
            lo, hi = lims[k]
            fmt = "{:.2f}" if k in ("voltage", "kinetics") else "{:.0f}"
            for yv, raw, va in ((1.0, hi / sign, "bottom"), (0.0, lo / sign, "top")):
                ax.text(xi, yv + (0.035 if va == "bottom" else -0.035), fmt.format(raw),
                        ha="center", va=va, fontsize=FS_TICK - 1, color="#333",
                        bbox=dict(fc="white", ec="none", pad=1.0), zorder=6)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{lab}\n({unit})" for _, lab, unit, _ in AXES], fontsize=FS_TICK)
        ax.set_ylim(-0.22, 1.22); ax.set_yticks([])
        ax.set_xlim(-0.35, len(AXES) - 0.65)
        ax.set_title(f"{s}\nfront: {', '.join(front_ids)}", fontsize=FS_TITLE, pad=12)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="x", length=0, pad=14)
    ids = sorted(cmap)
    fig.legend(handles=_legend_handles(ids, cmap), loc="center left",
               bbox_to_anchor=(0.995, 0.5), frameon=False, fontsize=FS_LEG)
    fig.text(0.01, 0.96, "↑ better on every axis (display-scaled; raw range at axis ends; "
             "whiskers = σ on front members, ▲▼ = σ exceeds the plotted range)",
             fontsize=FS_TICK - 1, color="#444")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = OUT / "pareto_parallel_coords.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    return p


# ---------------------------------------------------------------- 2. scatter matrix
def scatter_matrix(data, cmap):
    paths = []
    n = len(AXES)
    for s, rows in data.items():
        fig, axs = plt.subplots(n - 1, n - 1, figsize=(5.8 * (n - 1), 5.4 * (n - 1)))
        for i in range(1, n):
            for j in range(n - 1):
                ax = axs[i - 1, j]
                if j >= i:
                    ax.axis("off"); continue
                kx, lx, ux, _ = AXES[j]; ky, ly, uy, _ = AXES[i]
                for r in sorted(rows, key=lambda r: r["front"]):
                    vx, vy = r["raw"][kx], r["raw"][ky]
                    if vx is None or vy is None:
                        continue
                    col = cmap[r["id"]]
                    ax.errorbar(vx, vy, xerr=r["sig"][kx] or None, yerr=r["sig"][ky] or None,
                                fmt="o", ms=14, mfc=col if r["front"] else "white", mec=col,
                                mew=2.4, ecolor=col, elinewidth=1.4, capsize=4,
                                alpha=1.0 if r["front"] else 0.9,
                                zorder=3 if r["front"] else 2)
                if SIGN[kx] < 0:
                    ax.invert_xaxis()
                if SIGN[ky] < 0:
                    ax.invert_yaxis()
                if i == n - 1:
                    ax.set_xlabel(f"{lx}\n({ux})", fontsize=FS_LABEL)
                if j == 0:
                    ax.set_ylabel(f"{ly}\n({uy})", fontsize=FS_LABEL)
                ax.grid(True, color=ps.C["grid"], lw=0.6, alpha=0.45)
                ax.set_axisbelow(True)
                _style_axes(ax)
        front_ids = sorted(r["id"] for r in rows if r["front"])
        leg_ax = axs[0, n - 2]
        leg_ax.legend(handles=_legend_handles(sorted(cmap), cmap), loc="center",
                      frameon=False, fontsize=FS_LEG, title=f"{s}\nfront: {', '.join(front_ids)}",
                      title_fontsize=FS_LEG + 1)
        fig.suptitle(f"Pairwise objectives — {s}  (axes oriented so up/right = better; "
                     f"bars = σ)", fontsize=FS_TITLE, y=1.0)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        p = OUT / f"pareto_scatter_matrix_{s}.png"
        fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
        paths.append(p)
    return paths


# ---------------------------------------------------------------- 3. dominance heatmap
def dominance_heatmap(data, cmap):
    scen = list(data)
    fig, axs = plt.subplots(len(scen), 1, figsize=(15.5, 7.6 * len(scen)),
                            gridspec_kw=dict(hspace=0.45))
    axs = np.atleast_1d(axs)
    im = None
    for ax, s in zip(axs, scen):
        rows = sorted(data[s], key=lambda r: (not r["front"], r["id"]))
        m = len(rows)
        rank = np.full((m, len(AXES)), np.nan)
        for j, (k, *_rest) in enumerate(AXES):
            vals = [(SIGN[k] * r["raw"][k], i) for i, r in enumerate(rows)
                    if r["raw"][k] is not None]
            for rk, (_, i) in enumerate(sorted(vals, key=lambda t: -t[0])):
                rank[i, j] = rk + 1
        im = ax.imshow(rank, cmap="RdYlGn_r", vmin=1, vmax=m, aspect="auto")
        for i, r in enumerate(rows):
            for j, (k, *_rest) in enumerate(AXES):
                v = r["raw"][k]
                txt = "—" if v is None else (f"{v:+.2f}" if k == "voltage" else
                                             f"{v:.2f}" if k == "kinetics" else f"{v:.0f}")
                if k == "kinetics" and v is not None and r["lam_qc"] not in ("", "ok"):
                    txt += "*"
                ax.text(j, i, txt, ha="center", va="center", fontsize=FS_CELL)
            lines = [", ".join(r["dominated_by"][q:q + 3])
                     for q in range(0, len(r["dominated_by"]), 3)]
            label = "on front" if r["front"] else "dominated by\n" + "\n".join(lines)
            ax.text(len(AXES) - 0.3, i, label, va="center", ha="left", fontsize=FS_TICK - 2,
                    linespacing=1.1, color="#1B7837" if r["front"] else "#444")
        ax.set_xticks(range(len(AXES)))
        ax.set_xticklabels([f"{lab}\n({unit})" for _, lab, unit, _ in AXES],
                           fontsize=FS_TICK - 1)
        ax.set_yticks(range(m))
        ax.set_yticklabels([r["id"] for r in rows], fontsize=FS_TICK)
        for tl, r in zip(ax.get_yticklabels(), rows):
            tl.set_color(cmap[r["id"]]); tl.set_fontweight("bold" if r["front"] else "normal")
        ax.set_title(s, fontsize=FS_TITLE, pad=12)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
    cb = fig.colorbar(im, ax=axs.tolist(), orientation="horizontal", fraction=0.025,
                      pad=0.07, aspect=45, shrink=0.6)
    cb.set_label("rank in pool (1 = best)", fontsize=FS_TICK)
    cb.ax.tick_params(labelsize=FS_TICK - 2)
    fig.suptitle("Why each candidate is on / off the front\n(cells: raw value; colour = "
                 "per-axis rank; bold = front; * = λ QC-partial)", fontsize=FS_TITLE - 2,
                 y=0.97)
    p = OUT / "pareto_dominance_heatmap.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    return p


def main():
    ps.apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    data = load()
    cmap = _palette(data)
    for p in [parallel_coords(data, cmap), *scatter_matrix(data, cmap),
              dominance_heatmap(data, cmap)]:
        print(f"wrote {p}")


if __name__ == "__main__":
    main()
