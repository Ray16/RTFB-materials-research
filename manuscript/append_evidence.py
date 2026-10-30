#!/usr/bin/env python
"""Freeze pre-reviewed additions; updating the current pointer remains explicit."""
import argparse
import json
from pathlib import Path
import subprocess

from redox.validation.evidence_epoch import append_epoch

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--name',required=True)
    parser.add_argument('--extras',type=Path,required=True);parser.add_argument('--scope',required=True)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    print(json.dumps(append_epoch(root,args.name,json.loads(args.extras.read_text()),args.scope,commit),indent=2))
