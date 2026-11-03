# Scope of the closest experimental viologen comparison

The [identity and overlap audit](review.json) checks the soluble mediator in
[Kim et al. (2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9310624/) against the
paper's local graft model. The article reports first/second midpoint potentials
of -0.785/-1.204 V versus Fc+/0 in MeCN with 0.50 M tetrabutylammonium PF6.
The mediator is the benzyl/methyl viologen bis-PF6 salt. These measurements
provide useful analogue context and remain valid retrospective benchmark data.

They do not supply a new exact-graft validation. The measured graph is
C18H18N2(2+); our 4-methylbenzyl local model is C19H20N2(2+). Canonical graph
comparison confirms the additional para-methyl substituent. The computed
continuum/separated-ion setting also differs from the experimental electrolyte
and explicit salt. Neither a small structural change nor agreement of one wave
establishes equivalence of those conditions.

Both measured endpoints already occur in the frozen 39-endpoint benchmark.
For that same measured parent, the retained calculated potentials are
-0.594615/-1.207346 V, giving signed errors of +0.190385/-0.003346 V. These are
existing benchmark errors, not errors for the distinct capped graft. The audit
therefore leaves `new_exact_graft_error_V` blank and adds zero independent
held-out endpoints. No fitted correction or benchmark result is changed.

An additional quantitative exact-graft test needs the matching structure and
declared measurement conditions, with independence from the existing evaluation
panel stated explicitly. Until then, retain the article's distinction between
route precedent, analogue evidence and exact-graft predictive accuracy. The
final editorial PDFs are unchanged. Run `audit.py` through the project Python
environment to inspect the identity/overlap reconstruction; it runs no QM.
