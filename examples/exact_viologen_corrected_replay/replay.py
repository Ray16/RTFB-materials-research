#!/usr/bin/env python3
"""Run the self-contained corrected-cohort numerical cache."""
from pathlib import Path
import runpy
import sys

if __name__=='__main__':
    target=Path(__file__).resolve().parent/'data/replay.py'
    sys.argv[0]=str(target)
    runpy.run_path(str(target),run_name='__main__')
