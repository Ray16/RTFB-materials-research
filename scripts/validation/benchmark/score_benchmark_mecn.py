#!/usr/bin/env python
"""Compatibility export; implementation is redox.reports.benchmark."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from redox.reports.benchmark import *  # noqa: F403

if __name__ == "__main__":
    main()
