#!/usr/bin/env python
"""Recompute saved native-Hessian partitions with explicit source/xTB relocation."""
import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xtb-prefix', type=Path, required=True)
    parser.add_argument('--original-source-root', type=Path, required=True)
    parser.add_argument('--only', required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.threads < 1:
        parser.error('Positive thread count required')
    if args.out.exists():
        raise FileExistsError('Preserve previous runs')
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        os.environ[key] = str(args.threads)
    from redox.core.io import atomic_write_json
    from redox.qm.xtb_runtime import qualified_xtb_prefix
    from redox.qm.relocated_native import relocated_native_reader
    from redox.qm import projected_thermal
    with qualified_xtb_prefix(args.xtb_prefix) as runtime:
        with relocated_native_reader(args.original_source_root) as relocation:
            projected_thermal.main(['--only', args.only, '--threads', str(args.threads), '--out', str(args.out)])
            report = json.loads((args.out / 'run.json').read_text())
            if not report['checked'] or report['failed'] or report['skipped']:
                raise RuntimeError('Every declared relocated native state must be checked')
            atomic_write_json(args.out / 'runtime_binding.json', dict(runtime=runtime, relocation=relocation,
                              native_records_modified=False, quantum_calculations_run=False))


if __name__ == '__main__':
    main()
