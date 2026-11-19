# Figure style contract — apply to EVERY manuscript figure

One source of truth so the whole figure set reads as one system. All figures use
`scripts/plotting/plot_style.py` `apply_style()` (18 pt floor, one panel per file,
300 dpi, CVD-safe). Refine to ALL of the following; re-render; then READ the PNG
back and confirm no text < 18 pt and no label/legend/marker collisions.

## 1. Abbreviations (UNIFY — this is the #1 cross-figure bug)
- **DMPZ** = 5,10-dihydrophenazine. ALWAYS "DMPZ". Never "DHP". (Matches the
  `dmpz` candidate IDs.)
- **PTZ** = phenothiazine. **AQ/NQ/BQ** = anthra-/naphtho-/benzoquinone.
  **PMDI/NDI** = the imides. Define any abbreviation once per caption.
- First use in each caption spells the name out, then the abbreviation.

## 2. Family → colour (IDENTICAL across all figures; from plot_style palette)
Assign by entity, never by rank; keep the SAME hue for a family in every figure:
- quinone (n-type)      → blue   `#0072B2`
- imide (n-type)        → green  `#009E73`
- viologen / pyridine-multi-e (n-type) → reddish-purple `#CC79A7`
- phenazine / DMPZ (p-type) → orange `#E69F00`
- phenothiazine / PTZ (p-type) → purple `#762A83`
- nitroxide (p-type, validation only) → yellow/orange `#E69F00`→ use `#F0E442` to
  stay distinct from DMPZ orange on the parity plot
Legend always present for ≥2 series; ≤4 series also direct-labelled.

## 3. Axis labels & units (verbatim, every figure)
- x potential: `Redox potential, E vs Fc/Fc⁺ (V)`
- y capacity:  `Usable capacity (mAh g⁻¹, TBAPF₆/MeCN)`
- Reference always written `Fc/Fc⁺`. Potentials: 2 decimals, explicit `+` sign
  for positive (e.g. `+0.47 V`, `−1.54 V` with a true minus `−`).

## 4. Validation benchmark numbers (never drift)
MAE 0.247 V, RMSE 0.279 V, bias +0.19 V, n = 39 (30 molecules, 7 studies).
Any ±MAE band/whisker labelled `±0.25 V (MAE 0.247)` so the rounded band and the
exact stat don't look contradictory. Never 0.226 / 44.

## 5. Scientific-integrity rules (from AGENTS.md)
- Real chemical names. NO weighted-sum / scalar score — potential and capacity
  stay separate axes; anolyte and catholyte pools stay separate figures.
- p-type = catholyte (oxidation); n-type = anolyte (reduction) — label consistently.
- Honest status: potentials/capacity validated; λ_i qualitative; cell voltages are
  COMPUTED PROPOSALS (E_cell = E_cat − E_an), not demonstrated performance;
  capped graft ≠ soluble parent ≠ oligomer. Flag provisional (live-campaign) points.

## 6. Output (every figure)
- **PNG ONLY, 300 dpi.** Do NOT write `.pdf` (and delete any existing `.pdf` in
  `manuscript/figures/`; generators must not call savefig for pdf).
- `manuscript/figures/<figN>_<name>.png`.
- `manuscript/figures/<figN>_<name>.caption.md` (what · how · honest caveat).
- Generator in `scripts/plotting/manuscript/<figN>_<name>.py`, data read live
  (nothing hardcoded that could drift), frozen-91 pool enforced where relevant.

## 7. Text size & overlaps (HARD requirements — user-mandated)
- **Every piece of text is 18 pt** (titles, axis labels, ticks, legend, annotations,
  in-figure notes, numbers) — 18 pt is the floor; nothing smaller, anywhere.
- **ZERO overlapping elements.** No overlap between: text↔text, text↔marker,
  legend↔data, annotation↔annotation, label↔axis. Use leader lines, relocation,
  padding, or a larger figure canvas as needed. After rendering, READ the PNG and
  visually confirm there is not a single collision before reporting done.
