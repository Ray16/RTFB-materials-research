"""Nuclear partitions and thermal postprocessing.

``xtb_diagnostic`` is the GFN2/ALPB diagnostic model. ``native`` is solution
DFT Hessian thermochemistry. ``projected`` recomputes partitions from admitted
Hessians. The certified qRRHO implementation remains ``redox.qm.qrrho`` so
that its certificate hash and sibling ``xtb_kernel.py`` stay valid. Import
that module directly; this package does not load it at import time.
"""
