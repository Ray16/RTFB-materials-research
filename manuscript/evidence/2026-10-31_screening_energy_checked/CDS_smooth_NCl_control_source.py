#!/usr/bin/env python
"""Repeat full model controls for pyramidal ammonia and chloromethyl radical."""
import inspect
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    source=ROOT/'scripts/validation/qualify_smooth_quantum.py'
    module=api.load_module('unchanged_smooth_full_quantum_validation',source)
    body=inspect.getsource(module.main)
    replacements={
        '4\\nSinglet SMD model control\\nC 0 0 0\\nO 1.21 0 0\\nH -.61 .94 0\\nH -.61 -.94 0\\n':
            '4\\nPyramidal ammonia SMD model control\\nN 0 0 .25\\nH .94 0 -.1\\nH -.47 .814 -.1\\nH -.47 -.814 -.1\\n',
        '4\\nDoublet SMD model control\\nC 0 0 0\\nH 1.08 0 0\\nH -.54 .935307436 0\\nH -.54 -.935307436 0\\n':
            '4\\nChloromethyl radical SMD model control\\nC 0 0 0\\nCl 1.75 0 0\\nH -.54 .935307436 0\\nH -.54 -.935307436 0\\n'}
    for old,new in replacements.items():
        if body.count(old)!=1:raise ValueError('Full four-atom control input template changed')
        body=body.replace(old,new,1)
    namespace=dict(vars(module),__file__=str(Path(__file__).resolve()));local={}
    exec(compile(body,str(Path(__file__).resolve()),'exec'),namespace,local)
    try:
        local['main']()
    finally:
        out=Path(sys.argv[sys.argv.index('--out')+1]);path=out/'status.json'
        if path.exists():
            row=json.loads(path.read_text());row['source_sha256'][str(source)]=api.digest(source)
            row['chemical_control']='ammonia' if int(sys.argv[sys.argv.index('--mult')+1])==1 else 'chloromethyl_radical'
            atomic_write_json(path,row)


if __name__=='__main__':main()
