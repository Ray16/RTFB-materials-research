# Separate CPU PCM integral-memory budget

The [first large CPU reference point](../pmdi_split_execution_2026-10-24/electronic/status.json)
stopped at its declared 300 GiB RSS limit, during `PCM._get_v` three-center
integral construction. Its owner has ended and no energy is admitted. PySCF
already tiles these contractions, but inherits a large memory setting for PCM;
successive integral buffers can overlap in lifetime. The separate molecule/DF
and PCM budgets therefore need explicit control.

The [adapter](../../../scripts/pipeline/cpu_pcm_memory.py) caps the PCM object's
existing `max_memory` setting and records actual integral-tile dimensions. The
installed integral, surface, SCF and contraction code remains unchanged. The
new profile uses the numeric PCM memory setting 8192 in PySCF's MB convention;
the first pilot schema retains a `_MiB` suffix on this API value. Host RSS bounds
are measured separately in bytes/GiB. The solver retains a 400-point minimum
tile, so the PCM setting is not a strict process-RSS ceiling.

[Two spin controls](controls/status.json) exercise both potential and matrix
contractions with multiple density/charge sets and forced 400-point tiles.
They use the exact saved geometries and SCF densities of the prior checked
CPU/GPU controls. Both potential differences are zero; maximum matrix
differences are 2.43e-16 and 2.08e-16, below the predeclared 1e-10 bound.
The actual four tiles per contraction are recorded. These tests establish
contraction equivalence, not large-case resource suitability or native G.

The [fresh full-size point](large_point/plan.json) keeps the original geometry,
ωB97M-V/def2-TZVPD gas/SMD protocol, 1e-5 eV GPU-reference tolerance and resource
limits. Its [actual status](large_point/electronic/status.json) remains separate
from preparation and control success. A successful point and completed budget
trace are required before the normal geometry-to-energy handoff can use it.
No original failed result, model parameter or physical threshold is replaced.

Source and numerical evidence are retained beside this guide. The
[qualification program](../../../scripts/validation/qualify_cpu_pcm_memory.py)
also preserves unmodified method extracts and cached matrix arrays for replay.
Large temporary integrals use task-owned local scratch; input, source, SCF and
resource records are checkpointed on shared storage.
