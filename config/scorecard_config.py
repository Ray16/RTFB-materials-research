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
        (see redox.validate_stability) -> sigma ~0.15 eV [provisional; update from validation]
  - capacity (n, MW, specific): EXACT -> sigma 0
  - solubility proxy / SA: RELATIVE/heuristic only -> rank-only, no absolute sigma

Trust levels:
  exact      : bookkeeping, no model error (capacity)
  validated  : benchmarked against experiment/independent data (redox ranking, lambda, disprop)
  proxy      : rank-only, not an absolute prediction (solubility, SA)
  flagged    : known-unreliable, annotation only (dimerization absolute)
"""

# redox potential sigma (V) by couple type
SIGMA_REDOX_V = dict(cation=0.44, anion=0.53, multi=0.80)
SIGMA_LAMBDA_EV = 0.10
SIGMA_DISP_EV = 0.17          # measured: redox.validate_stability MAE 16 kJ/mol (n=3, rho=1.0)
SIGMA_CAPACITY = 0.0          # exact

TRUST = dict(
    redox_potential="validated",     # ranking within class; absolute carries the sigma above
    reorganization="validated",
    disproportionation="validated",
    capacity="exact",
    solubility="proxy",
    synthetic_accessibility="proxy",
    dimerization="flagged",
)


def sigma_redox(q_ox, q_red):
    """Redox-potential sigma (V) for a couple, from the OROP by-charge-class benchmark."""
    if abs(int(q_ox)) >= 2 or abs(int(q_red)) >= 2:
        return SIGMA_REDOX_V["multi"]
    # a couple that produces/consumes an anion (min charge < 0) is a reduction couple
    return SIGMA_REDOX_V["anion"] if min(int(q_ox), int(q_red)) < 0 else SIGMA_REDOX_V["cation"]
