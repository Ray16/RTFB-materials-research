"""Import path retained for the hash-locked occupied-sampling figure script.

The implementation is ``redox.execution.resources``.
"""
from redox.execution.resources import *  # noqa: F401,F403
from redox.execution import resources as _implementation
import sys
sys.modules[__name__] = _implementation
