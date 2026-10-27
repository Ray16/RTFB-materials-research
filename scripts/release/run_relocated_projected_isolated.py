#!/usr/bin/env python
"""Run the relocated projected CLI with original-checkout Python opens blocked."""
import argparse
import os
from pathlib import Path
import runpy
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('project','forbid-root','xtb-prefix','original-source-root','site-packages','out'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--only',required=True)
    args=parser.parse_args()
    if not sys.flags.no_site:
        parser.error('Run with -S to avoid editable installation hooks')
    root=args.project.resolve();forbidden=args.forbid_root.resolve()
    if root.is_relative_to(forbidden) or args.xtb_prefix.resolve().is_relative_to(forbidden):
        parser.error('Relocated project and xTB prefix must be outside original checkout')
    sys.path.append(str(args.site_packages.resolve()));sys.path.insert(0,str(root/'src'))
    def audit(event,values):
        if event=='open' and values and isinstance(values[0],(str,bytes,os.PathLike)):
            path=Path(os.fsdecode(values[0])).resolve()
            if path.is_relative_to(forbidden):
                raise PermissionError('Original checkout access blocked: '+str(path))
    sys.addaudithook(audit)
    entry=root/'scripts/pipeline/run_relocated_projected.py'
    sys.argv=[str(entry),'--xtb-prefix',str(args.xtb_prefix),'--original-source-root',str(args.original_source_root),
              '--only',args.only,'--threads','1','--out',str(args.out)]
    runpy.run_path(str(entry),run_name='__main__')


if __name__=='__main__':main()
