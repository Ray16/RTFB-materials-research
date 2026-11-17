"""Point child processes at this checkout rather than another installed copy."""
import os
from pathlib import Path

_SRC = str(Path(__file__).resolve().parents[1] / "src")
_CURRENT = os.environ.get("PYTHONPATH")
os.environ["PYTHONPATH"] = _SRC if not _CURRENT else _SRC + os.pathsep + _CURRENT
