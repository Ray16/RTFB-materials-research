# Architecture

Dependencies run from `core` to the scientific packages, then to orchestration
and reporting. Backends do not import campaign orchestration. Optional
quantum libraries load on the paths that call them.

## Layers

- **structures** builds exact precursor/product grafts, capped models, soluble
  parents and open oligomers. It checks graph identity, charge, spin and
  stereochemistry. It does not score a route or a potential.
- **backends** holds DFT/SMD, UMA selection, solvent response and checkpoint
  IO. Hessian recovery stays four explicit modules: `hessian_memory` records
  an execution budget, `hessian_resume` reconstructs an additive CDS term,
  `hessian_restart` reuses an orbital frame, and `saved_hessian` reuses one
  complete matrix. `response_resume` may reuse a response only after a fresh
  equation check. `response_diagnostic` is instrumentation and does not
  qualify a Hessian.
- **thermochemistry** separates the GFN2/ALPB diagnostic (`xtb_diagnostic`),
  solution DFT Hessian thermochemistry (`native`), its minimum-repair
  promotion (`native_repair`), rigid-mode projection (`vibrations`) and
  projected partitions (`projected`). These are different thermal definitions.
- **qm** keeps the certified qRRHO kernel. `qrrho.py` and the sibling
  `xtb_kernel.py` are byte-identical to the source tree. The projected-thermal
  certificate admits only one source transformation, the xTB library path, and
  locates the kernel worker beside `qrrho.py`. `xtb_runtime.py` stays
  importable as `redox.qm.xtb_runtime` because that import is inside the
  certified bytes.
- **properties** computes potentials, reorganization, grafting thermodynamics,
  ensembles and capacity. Half-reactions and sums stay separate. There is no
  scalarized figure of merit.
- **screening** registers exact precursor/product pairs, applies
  proceed/hold/exclude route gates, and keeps Pareto selection over raw axes.
- **quality** reads numerical, physical, source and reference evidence. A
  finished checkpoint is not a qualified free energy.
- **execution** discovers GPUs through `gpu_reserve`, caps threads, and moves
  checkpoint ownership. It does not change equations.
- **workflows** are campaign and recovery policies. Shared stage execution is
  `workflows.stages`. Resource acquisition is not redefined in each campaign.
- **reporting** rebuilds tables from admitted records.

## Data flow

Exact precursor/product definitions and evidence, then route
proceed/hold/exclude, then state and oligomer structures, then
protocol-addressed UMA/DFT geometry and energy records, then source-bound
thermal and QC admission, then matched-reference potentials, reorganization,
reaction free energies and finite-occupancy quantities. Pending values stay
absent.

## Certified bytes

Moving a producer changes its hash. Historical certificates compare those
hashes. This tree therefore does not rewrite:

- `src/redox/qm/qrrho.py`
- `src/redox/qm/xtb_kernel.py`
- `src/redox/quality/sampling_snapshot.py`, whose hash is the frozen
  occupied-sampling producer
- `scripts/plotting/plot_occupied_sampling.py`, whose hash is the frozen figure
  producer
- hash-locked files under `examples/`

`redox.workflows.resources` and `redox.validation.sampling_receipts` /
`sampling_native_receipts` are compatibility imports for those frozen
producers. Their implementations live in `execution.resources` and `quality`.
