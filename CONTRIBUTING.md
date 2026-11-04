# Contributing

Use the canonical project checkout and preserve unrelated edits. Inspect
`python -m redox workflows` to distinguish production, analysis, control and
recovery commands. Run commands from the repository root with the appropriate
environment. Python configuration files execute code and must be trusted.

## Small changes

Explain the concrete problem, expected behavior, changed interfaces and checks.
Prefer shared calculation/validation logic over implementations that can diverge.
Keep changes focused and preserve existing interfaces where practical.

A lightweight CPU check in an environment with pytest is:

```bash
PYTHONPATH=src python -m pytest -q tests/test_cached_thermochemistry.py tests/test_workflow_catalog.py tests/test_config_loading.py
python3 -I -S examples/cached_thermochemistry/reproduce.py
```

The complete suite is `python -m pytest -q` in the scientific environment. Report
exactly what ran, including skipped/failed tests. Hosted CI covers only the small
subset above. It does not validate CUDA execution, Hessians or experimental accuracy.

## Scientific changes

Preserve equations, units, charge/spin, solvent, reference, thermal definitions,
random seeds and QC thresholds unless the change explicitly requires otherwise.
Record before/after results with justified tolerances. New numerical methods,
architectures or runtime versions need an independently scoped validation gate.
Do not convert missing components to zero or promote a cached replay into new
physical qualification. Preserve source-frozen records and historical results.

Keep raw data/checkpoints under the canonical project's existing data/output
locations. Include hashes, configuration, software/hardware context and clear
prepared/pending/qualified status in small reviewable records. Do not commit
credentials, private data, downloaded weights or bulk calculation artifacts.

## Reporting problems

Include the command, exact commit, environment versions, relevant configuration,
input identity, expected/observed behavior and a minimal shareable example.
Remove tokens and personal information from logs. For calculation failures,
preserve checkpoints and the original error before attempting recovery.

The existing MIT license applies to repository software. Contributors must
identify separate restrictions on third-party data, source papers and weights.
Do not invent author metadata or a publication citation.
