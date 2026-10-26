#!/usr/bin/env python
"""Run the relocated xTB kernel probe with original-project Python opens blocked.

Invoke with Python -S to avoid editable-install .pth files. Trace child-process
opens separately if system-call-level evidence is needed.
"""
import argparse
import os
from pathlib import Path
import runpy
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--forbid-root', type=Path, required=True)
    parser.add_argument('--xtb-prefix', type=Path, required=True)
    parser.add_argument('--site-packages', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if not sys.flags.no_site:
        parser.error('Run Python with -S to disable editable-install .pth files')
    project, forbidden = args.project.resolve(), args.forbid_root.resolve()
    if project.is_relative_to(forbidden) or args.xtb_prefix.resolve().is_relative_to(forbidden):
        parser.error('Relocated source and xTB prefix must be outside forbidden root')
    sys.path.append(str(args.site_packages.resolve()))
    sys.path.insert(0, str(project / 'src'))

    def audit(event, values):
        if event == 'open' and values and isinstance(values[0], (str, bytes, os.PathLike)):
            path = Path(os.fsdecode(values[0])).resolve()
            if path.is_relative_to(forbidden):
                raise PermissionError(f'Original project access blocked: {path}')

    sys.addaudithook(audit)
    entrypoint = project / 'scripts/pipeline/run_portable_thermal.py'
    sys.argv = [str(entrypoint), '--mode', 'probe', '--xtb-prefix', str(args.xtb_prefix),
                '--out', str(args.out)]
    runpy.run_path(str(entrypoint), run_name='__main__')


if __name__ == '__main__':
    main()
