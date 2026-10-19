"""Discoverable entry points for researchers and agents: ``python -m redox --help``."""
from __future__ import annotations

import argparse
import runpy
import sys


COMMANDS = {
    "status": ("redox.validation.readiness", (), "inspect all declared states and QC"),
    "energy-audit": ("redox.qm.sp", ("--audit",), "inspect active single-point coverage; no QM"),
    "review": ("redox.workflows.tables", (), "rebuild and audit tables in a fresh directory; no QM"),
    "thermal-postprocess": ("redox.qm.projected_thermal", (), "recompute corrected partitions from admitted Hessians; CPU only"),
    "build": ("redox.build.candidates", (), "build candidate structures; --dry-run validates only"),
    "expand": ("redox.screening.expansion", (), "prepare matched-comparison hypotheses; --register appends definitions"),
    "search": ("redox.screening.family_search", (), "enumerate family/substituent/linker search space; no QM or registration"),
    "source-audit": ("redox.screening.source_audit", (), "audit downloaded source quality; no label admission or QM"),
    "resources": ("redox.workflows.resources", (), "discover free GPUs through the system policy"),
    "launch": ("redox.workflows.launch", (), "plan DFT/UMA shards; --execute submits"),
    "run-gpu": ("redox.workflows.gpu", (), "preflight and run one reserved full-QM command"),
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
