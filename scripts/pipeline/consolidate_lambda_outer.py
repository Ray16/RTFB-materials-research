#!/usr/bin/env python
"""Compatibility export; implementation is redox.reports.outer_sphere."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from redox.reports.outer_sphere import *  # noqa: F403

if __name__ == "__main__":
    main()
