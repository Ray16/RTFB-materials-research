# Standalone route-gate identity guard

The standalone prospective route screen previously matched evidence by ID and
could return `proceed` for a changed fragment without independently binding the
review to that fragment. The streamed-route wrapper and all-axis decision
collector already performed separate identity checks. The standalone gate now
reconstructs the capped-site product, precursor and mapped new bond, and permits
`proceed` only when they equal the reviewed binding. Declared precursor and
product graphs must agree as well. Placeholder parity-valid states used for
identity construction are not redox-state predictions. This guard does not
validate the scientific claims or sources in the review.

Checks: the direct gate accepts a matching mapped route, rejects a changed bond,
missing binding and conflicting product; the existing TEMPO exact-support route
still proceeds with an exact match. The default 75-candidate screen remains at
75 holds. Focused cached thermochemistry, topology and direct-gate checks passed
(10 tests); streamed-route checks passed (9 tests). The earlier broader focused
run passed 26 tests, including archived polymer replay; it was started before
this guard was edited, so it is not a post-change regression for the new guard.
No quantum calculations were run. The previous audit's gate hash and output are
historical snapshots, not a claim that the gate source is unchanged.
