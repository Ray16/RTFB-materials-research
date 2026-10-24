# V100 UMA execution and model-version diagnosis

**V100 support is retained.** The production profile remains FairChem 2.21.0,
Torch 2.8.0+cu128 and the fingerprinted UMA-S-1.2.1 checkpoint. Its ten-case,
forty-frame [legacy panel](../uma_nitro_2026-10-23/legacy_points.json) runs on a
reserved V100; all six fingerprinted stored energies reproduce within
6.10e-7 eV. These are gas-phase UMA tests. They do not certify an arbitrary
polymer DFT/Hessian calculation or replace its physical validation.

The newer FairChem commit `645af7d03e91138ef07e51de387316c9e2e152dc`
also executes all ten states on V100 with Torch 2.13.0+cu126. The
[architecture probe](architecture_probe.json) confirms `sm_70`, actual CUDA
arithmetic and the V100 device. The [explicit FP32 worker](../../../scripts/validation/qualify_uma_v100.py)
sets `nitro=False`, `tf32=False`, `compile=False`, `merge_mole=False`,
`activation_checkpointing=True`, graph version 2 and `execution_mode='general'`.
All 290 post-evaluation mode checks pass. Both reservations ended normally,
without foreign-process interruption or modifications to installed libraries.

Execution alone does not establish reproduction: the [unaltered current-model
comparison](reference_comparison.json) passes only 3/10 cases under the frozen
energy **and** force bounds. Differences reach 0.0274357 eV and
0.138183 eV/Å. The [cached points](reference_points.json) retain every frame,
not just averages or selected neutral states. Original tolerances are unchanged.

The [source comparison](source_snapshot.json) identifies a concrete behavior
change. The legacy loader uses `model_version=1.0`, including a zero-initialized
accumulator in the composition-embedding mean; the newer 1.2.1 implementation
excludes that accumulator. These give sums divided by N+1 and N, respectively.
An [isolated diagnostic](../../../scripts/validation/diagnose_uma_composition.py)
changes only the in-memory averaging branch to the historical rule. All
[10 cases/40 frames then pass](composition_comparison.json), with maximum errors
4.59e-6 eV and 3.34e-6 eV/Å. The weights, atomic references, explicit charge/spin
and coordinates remain identical. This localizes the tested drift to composition
averaging. The [diagnostic points](composition_points.json) are explicitly marked
ineligible for production; its temporary model-ID override must not be adopted
as an official current-model configuration.

The newer CUDA13 Nitro prefix remains separate. Its Torch wheel does not include
`sm_70`; disabling Nitro cannot make that wheel support V100. The
[official API documentation](https://facebookresearch.github.io/fairchem/ase-calculator/)
also limits actual Nitro to compute capability 8.x–12.x. V100 therefore uses a
CUDA12 standard inference profile. Never label such execution as Nitro or report
a Nitro speedup from a V100 comparison.

The [summary](summary.json) records the admission limits and remaining migration
work. UMA-derived gas diagnostics, spin scans and warm starts must retain their
loader behavior and separate caches. Current FP32 versus actual Nitro redox
validation remains separate. No DFT energy, MeCN free energy or candidate decision
was replaced by these gas-phase tests.

The [CPU comparison script](../../../scripts/validation/review_uma_points.py)
recomputes every energy and force check from stored arrays and verifies panel
identities. The [tests](parity_tests.xml) include incorrect forces with matching
energies, missing geometries, altered spin and altered atomic references. On the
original evidence filesystem, reproduce the comparison with:

```bash
PYTHONPATH=src python scripts/validation/review_uma_points.py \
  --plan results/validation/uma_nitro_2026-10-23/panel/plan.json \
  --reference results/validation/uma_nitro_2026-10-23/legacy_points.json \
  --candidate results/validation/uma_v100_2026-10-24/reference_points.json \
  --out /tmp/uma_v100_reference_review.json
```

Raw evaluations remain in `calcs/uma_nitro_validation/2026-10-24_v100_reference_02/`
and `calcs/uma_nitro_validation/2026-10-24_composition_ablation_01/`.
Launch manifests, finite residency guards and the first editable-registry lookup
failure are preserved here. That first failure occurred before prediction; the
local identity reader now supports both wheel and editable layouts without
modifying their registries. The CUDA12 compatibility prefix was inspected and used
read-only; it is not the canonical shared Nitro environment.
