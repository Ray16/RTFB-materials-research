"""Certified thermodynamic kernel.

``qrrho.py`` and its sibling ``xtb_kernel.py`` stay here because the projected
thermal certificate admits only one source transformation, and the kernel worker
is located beside ``qrrho.py``. ``xtb_runtime.py`` stays importable at this path
because that import line is part of the certified ``qrrho.py`` bytes.

DFT, UMA, Hessians and solvent models live in ``redox.backends``. Nuclear
postprocessing lives in ``redox.thermochemistry``.
"""
