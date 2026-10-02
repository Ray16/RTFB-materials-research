"""Per-axis uncertainty (sigma) and trust levels for the candidate scorecard.

These are the numbers the selection engine needs to do sigma-aware domination and to know
which axes are trustworthy objectives vs proxies/filters. Values are grounded in our own
benchmarks (not guessed):

  - redox potential sigma: OROP experimental MAE, BY CHARGE CLASS
        cations (|q|<=1, oxidation side)   ~0.44 V
        anions  (|q|<=1, reduction side, reversible) ~0.53 V
        multiply-charged (|q|>=2)          ~0.80 V   (few points; conservative)
  - reorganization energy lambda: D3TaLES cross-check MAD ~0.094 -> sigma 0.10 eV
  - disproportionation dG: cancellation-limited; validated vs experimental wave spacing
        (see redox.validation.stability) -> sigma ~0.15 eV [provisional; update from validation]
  - capacity (n, MW, specific): EXACT -> sigma 0
  - solubility proxy / SA: RELATIVE/heuristic only -> rank-only, no absolute sigma

Trust levels:
  exact      : bookkeeping, no model error (capacity)
  validated  : benchmarked against experiment/independent data (redox ranking, lambda, disprop)
  proxy      : rank-only, not an absolute prediction (solubility, SA)
  flagged    : known-unreliable, annotation only (dimerization absolute)
"""

# redox potential sigma (V) by couple type — CONSERVATIVE OROP by-charge-class PRIOR, used
# only as a last-resort fallback for a family with NO validation data at all. The sigma the
# scorecard actually reports is DERIVED AT RUN TIME from our own validation residuals
# (see redox.screening.scorecard.derive_family_sigma): per-family RMSE where we have >=2 anchors, else
# the pooled validation RMSE (our overall demonstrated accuracy). Nothing is hard-coded, so
# the uncertainty tightens automatically as config/validation.py grows.
# NOTE: every SIGMA below is a FALLBACK only. The scorecard computes each uncertainty at run
# time from its benchmark file (redox.screening.scorecard.derive_uncertainties): redox from the OROP
# benchmark + validation residuals, lambda from the D3TaLES comparison, disproportionation
# from stability_validation. These constants are used ONLY if that benchmark file is absent,
# so no reported uncertainty is a hand-set magic number.
SIGMA_REDOX_V = dict(cation=0.44, anion=0.53, multi=0.80)   # fallback; else OROP-derived
SIGMA_LAMBDA_EV = 0.10        # fallback; else RMSE(our lambda - D3TaLES) on the 0-2 eV subset
SIGMA_DISP_EV = 0.17          # fallback; else RMSE of stability_validation residuals
SIGMA_CAPACITY = 0.0          # exact (bookkeeping, no model error)

# Outer-sphere lambda_o (molecular-cavity nonequilibrium PCM, 1-body, MeCN; results/lambda_outer.csv
# from scripts/pipeline/compute_lambda_outer.py). Continuum-vs-explicit gap is the dominant error: up to
# 0.42 eV for naphthalene/THF (Ambrosio et al. JPCL 2025, 10.1021/acs.jpclett.5c01328), expected
# smaller in polar MeCN. PRIOR (one benchmark; no MeCN explicit/experimental anchor yet), so
# lambda_o / lambda_total are RANKING-ONLY. See results/reorg_anchors/outer_sphere_validation.md.
SIGMA_LAMBDA_O_EV = 0.30

# Counter-ion-inclusive capacity, MAX-LOAD convention (symmetric for cationic and anionic states):
# repeat-unit mass + the HEAVIEST counter-ion load over all accessible states (q>0 -> q x PF6-,
# q<0 -> |q| x supporting cation). A "resting-state only" convention was rejected: it silently
# favours neutral-resting n-type molecules (their dianion cations would be ignored).
COUNTERION_ANION_SMILES = "F[P-](F)(F)(F)(F)F"   # PF6- (config/project.json "counterion")
SUPPORTING_CATIONS = {                            # supporting-salt cation per scenario
    "Li": "[Li+]",                                # LiPF6
    "TBA": "CCCC[N+](CCCC)(CCCC)CCCC",            # TBAPF6 (standard MeCN salt)
}
# Pareto capacity objective per salt scenario. The salt is not decided yet (user: report both,
# choose later), so redox.screening.pareto builds the front under EACH and reports where it changes.
CAPACITY_SCENARIOS = {
    "LiPF6": "capacity_maxload_Li_mAh_g",
    "TBAPF6": "capacity_maxload_TBA_mAh_g",
}

# Spin-state confidence threshold: flag a molecule when the smallest |spin-state gap| across
# its UMA states is near-degenerate (the DFT charge/energy could rest on the wrong
# multiplicity). This is a DEFINITIONAL cutoff (a modelling choice), not a benchmark number,
# so it stays here — expressed as its physical value, 5 kcal/mol, not a bare 0.217.
KCAL_PER_MOL_EV = 0.0433641   # eV per kcal/mol
SPIN_GAP_LOW_EV = 5.0 * KCAL_PER_MOL_EV

TRUST = dict(
    redox_potential="validated",     # ranking within class; absolute carries the sigma above
    reorganization="validated",
    disproportionation="validated",
    capacity="exact",
    reorganization_outer="proxy",    # continuum lambda_o: ranking only (sigma above)
    solubility="proxy",
    synthetic_accessibility="proxy",
    dimerization="flagged",
)


def sigma_redox(q_ox, q_red, family=None, family_sigma=None, pooled_sigma=None):
    """Redox-potential sigma (V) for a couple, in priority order:
      1. the per-FAMILY RMSE measured from our validation residuals (family_sigma), if present;
      2. the POOLED validation RMSE (pooled_sigma) — our overall demonstrated accuracy — for a
         family we have not directly validated (e.g. imide, until an imide anchor is added);
      3. only if we have NO validation at all, the conservative OROP by-charge-class prior.
    family_sigma / pooled_sigma are computed at run time (redox.screening.scorecard.derive_family_sigma),
    so nothing here is hard-coded and the uncertainty tightens as the validation set grows."""
    if family_sigma and family in family_sigma:
        return family_sigma[family]
    if pooled_sigma is not None:
        return pooled_sigma
    if abs(int(q_ox)) >= 2 or abs(int(q_red)) >= 2:
        return SIGMA_REDOX_V["multi"]
    # a couple that produces/consumes an anion (min charge < 0) is a reduction couple
    return SIGMA_REDOX_V["anion"] if min(int(q_ox), int(q_red)) < 0 else SIGMA_REDOX_V["cation"]
