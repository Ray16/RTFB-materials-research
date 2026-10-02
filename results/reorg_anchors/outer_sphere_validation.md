# Outer-sphere (solvent) reorganization energy — validation status

Goal: validate our λ_out (and hence λ_total = λ_in + λ_out) against literature, the same way
we validated λ_in (gas 4-point) to meV against Deng & Goddard / da Silva Filho / Delgado.

λ_in is validated (see oligoacene_* and perfluoroacene_* figures). This file tracks λ_out.

## Benchmark 1 — naphthalene in THF (OA, rigorous, extracted from full text)
Ambrosio, Landi, Loriso, Leo, Peluso, *J. Phys. Chem. Lett.* 2025 — "External Reorganization
Energy upon Charge Transfer Reactions in Mildly Polar Media: Naphthalene in THF."
DOI 10.1021/acs.jpclett.5c01328 (Europe PMC PMC12235641, open access).

Values reported (THF; naphthalene/naphthalene•− reduction):
- λ_int (inner, gas)              = 0.15 eV
- λ_ext(red), **continuum (PCM)** = **0.70 eV**   (1-body / half-reaction)
- λ_ext, self-exchange (2-body)   ≈ 1.1 eV
- **continuum vs explicit-solvent (thermodynamic integration): differ by up to 0.42 eV**
  (continuum neglects explicit solute–solvent interactions)

Our pipeline (Born single-sphere, SASA radius a=3.67 Å, THF eps_s=7.58 eps_op=1.98):
- **our λ_out (1-body) = 0.732 eV**  vs  paper PCM 0.70 eV  →  **agreement 0.03 eV**

### What this validates / tells us
1. **Our Born λ_out is a FAITHFUL CONTINUUM estimate** — it reproduces a molecular-cavity PCM
   value to ~0.03 eV. The SASA-radius choice (probe ~1.3 Å) is what makes it match; the bare
   vdW-volume radius overestimates (0.86 eV).
2. **Our Born value is the 1-body (electrochemical / self-exchange-at-contact) λ_out.** To compare
   to experimental SELF-EXCHANGE rates, use the 2-body form (~2x at infinite separation, reduced
   at contact — the paper's 0.70 -> 1.1 shows the contact reduction, not a clean 2x).
3. **The dominant error is continuum-vs-explicit (~0.42 eV here), NOT sphere-vs-molecular-cavity.**
   A fancier continuum (nonequilibrium molecular-cavity PCM) would ~reproduce the 0.70 eV, so it
   would NOT close this gap. Closing it requires explicit/cluster solvent. NB: THF is mildly polar
   (worst case for continuum); in MeCN (eps_s=37.5) the continuum-vs-explicit gap is likely smaller.

## Pending — MeCN experimental self-exchange λ (our families)
Needed to measure the pipeline's ACTUAL error in our solvent. Specific classic sources (paywalled;
digits to be retrieved):
- TEMPO•/TEMPO+ self-exchange, MeCN — Grampp & Rasmussen, PCCP 2002, 10.1039/b206313a
- Quinone/semiquinone self-exchange, MeCN — Grampp/Landgraf/Rasmussen 1999 (10.1039/a903394g);
  Williams 1969 (10.1080/00268976900100071)
- Phenothiazine/PTZ+ (IV-CT + self-exchange) — Sun, Rosokha, Kochi, JACS 2004, 10.1021/ja038746v
Compare to λ_total = λ_in(gas 4-point) + λ_out(Born, 2-body). Expect ~0.1-0.3 eV agreement
(continuum floor), NOT meV.

## Pitfalls being respected (per user)
1. No double-counting: λ_in is GAS-phase 4-point -> orthogonal to a vertical λ_out. Keep it gas.
2. 1-body vs 2-body: our Born = 1-body; use 2-body for self-exchange comparisons.
3. Optical dielectric: MeCN eps_op = n^2 = 1.344^2 = 1.806 (config/electrolyte.py); the Born
   formula uses eps_op explicitly. (PySCF PCM does NOT cleanly expose eps_op for ground-state
   nonequilibrium — a manual fast/slow partition would be required if we go the PCM route.)
