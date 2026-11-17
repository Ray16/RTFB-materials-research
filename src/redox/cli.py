"""Discoverable entry points for researchers and agents: ``python -m redox --help``."""
from __future__ import annotations

import argparse
import runpy
import sys


COMMANDS = {
    "workflows": ("redox.workflows.catalog", (), "show maintained workflows, calculation stages, controls and recovery tools; no execution"),
    "uma": ("redox.backends.uma_current", (), "gas-phase warm starts with explicit backend qualification"),
    "dft": ("redox.backends.uma_current", ("--stage", "dft"), "DFT/SMD with checked warm-start provenance"),
    "motifs": ("redox.workflows.motifs", (), "plan or execute registered candidates"),
    "uma-check": ("redox.quality.uma_environment", (), "inspect UMA configuration without loading CUDA"),
    "status": ("redox.quality.readiness", (), "inspect all declared states and QC"),
    "sp": ("redox.backends.sp", (), "active-protocol single points; adopts matching DFT energies (--backend gpu)"),
    "native-thermal": ("redox.thermochemistry.native", (), "solution DFT Hessian thermochemistry via gpu_reserve (--only, --select-policy)"),
    "seed-reference": ("redox.workflows.reference", (), "install + verify pinned Fc/benchmark calibration records (data/reference)"),
    "energy-audit": ("redox.backends.sp", ("--audit",), "inspect active single-point coverage; no QM"),
    "review": ("redox.reporting.tables", (), "rebuild and audit tables in a fresh directory; no QM"),
    "thermal-postprocess": ("redox.thermochemistry.projected", (), "recompute corrected partitions from admitted Hessians; CPU only"),
    "build": ("redox.structures.candidates", (), "build candidate structures; --dry-run validates only"),
    "build-validation": ("redox.structures.validation", (), "build validation cores incl. the ferrocene reference (--set benchmark)"),
    "expand": ("redox.screening.expansion", (), "prepare matched-comparison hypotheses; --register appends definitions"),
    "search": ("redox.screening.family_search", (), "enumerate family/substituent/linker search space; no QM or registration"),
    "source-audit": ("redox.screening.source_audit", (), "audit downloaded source quality; no label admission or QM"),
    "resources": ("redox.execution.resources", (), "discover free GPUs through the system policy"),
    "launch": ("redox.execution.launch_current", (), "plan current DFT/UMA shards; --execute submits"),
    "run-gpu": ("redox.execution.gpu", (), "preflight and run one reserved full-QM command"),
}


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description="Redox flow battery property pipeline")
    parser.add_argument("command", choices=COMMANDS)
    parser.epilog = "\n".join(f"{name:14s} {definition[2]}" for name, definition in COMMANDS.items())
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    if not args or args[0] in ("-h", "--help"):
        parser.print_help()
        return 0
    if args[0] not in COMMANDS:
        parser.error(f"unknown command {args[0]!r}")
    module, prefix, _ = COMMANDS[args[0]]
    previous = sys.argv
    try:
        sys.argv = [module, *prefix, *args[1:]]
        runpy.run_module(module, run_name="__main__")
    finally:
        sys.argv = previous
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
