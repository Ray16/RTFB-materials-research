"""Quantum backends.

DFT/SMD, UMA, solvent response and checkpoint IO live here. Four Hessian
strategies stay separate modules because they are different scientific
operations: ``hessian_memory`` (execution budget), ``hessian_resume`` (additive
CDS reconstruction), ``hessian_restart`` (orbital frame restart) and
``saved_hessian`` (complete-matrix reuse). Response reuse and response
diagnostics likewise stay separate.
"""
