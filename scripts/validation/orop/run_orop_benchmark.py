#!/usr/bin/env python
"""Compatibility entry point for the matched-reference OROP benchmark."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from redox.validation.orop import main


if __name__ == "__main__":
    raise SystemExit(main())
