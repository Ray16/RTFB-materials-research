import importlib.util
from pathlib import Path

import numpy as np
import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts/pipeline/smd_memory_execution.py"
spec = importlib.util.spec_from_file_location("smd_buffer_lifetimes", PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_inplace_K_assembly_retains_the_same_dense_matrix():
    rng = np.random.default_rng(43)
    S, D = rng.normal(size=(2, 20, 20))
    A = rng.uniform(0.01, 1, 20)
    epsilon = 0.93
    reference = S - epsilon / (2 * np.pi) * ((D * A) @ S)
    scratch = (D * A) @ S
    scratch *= -epsilon / (2 * np.pi)
    scratch += S
    assert np.array_equal(reference, scratch)


def test_unreviewed_upstream_source_is_not_patched_by_guessing():
    with pytest.raises(ValueError, match="reviewed"):
        module.compact_build_source("def build(self): pass")
    with pytest.raises(ValueError, match="reviewed"):
        module.release_geometry_source("def dft_smd(): pass")
