#!/usr/bin/env python
"""Compatibility entry point for ``python -m redox.workflows.tables``."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from redox.workflows.tables import main


if __name__ == "__main__":
    raise SystemExit(main())
