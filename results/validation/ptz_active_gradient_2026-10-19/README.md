# Active-level force diagnostic at the composite geometries

Three independent reserved H100 jobs reproduce selected neutral/radical-cation/
dication active energies at frozen r2SCAN-D4 geometries, using tighter SCF settings.
Full NLC and SMD grid-response gradients are checked against cold finite-energy
differences for the largest component at +/-0.0005 Bohr. All checks pass (absolute
errors below 2e-6 Ha/Bohr); maximum active gradients are about 0.015 Ha/Bohr.

This establishes nonstationarity on the active electronic surface, not a software
failure of the declared composite protocol. A lower-level native-Hessian local
quadratic estimate predicts energy changes of -0.080/-0.090/-0.109 eV, with a
0.1-A trust radius. Their estimated charging-gap effect is only -0.0084 eV versus
the observed 0.354-eV gap discrepancy. This is not an active-level minimization,
free-energy correction or rigorous bound. No proposed geometry, energy or thermal
correction is adopted and canonical inputs remain unchanged.

`plan.json`, node preflight, launch manifest and actual-start verification bind
source identities and resource use. Raw gradient arrays stay on shared storage
with hashes; small reports are versioned. The measured discrepancy remains open.

The [independent record review](independent_review.json) passes 63 checks,
including raw finite-energy quotients, gradient and source hashes, selected-energy
reproduction and the charging-gap identity. It does not repeat quantum jobs.
[Terminal verification](terminal_verification.json) confirms all three actual
workers and their controllers are absent; each reservation wrapper exited zero.
