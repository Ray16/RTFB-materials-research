"""Atomic replacement for records read by resumable workers and watchers.

Stage in the destination directory so replacement stays on the same filesystem.
Callers that merge shared state must also hold their domain's writer lock.
"""
from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import stat
import tempfile


@contextmanager
def atomic_text_writer(path: Path, *, newline=None):
    """Publish a complete UTF-8 file on success; retain the old file on failure."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    staged = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline=newline, dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as stream:
            staged = Path(stream.name)
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        # Preserve existing collaboration permissions; new records are readable.
        mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else 0o644
        staged.chmod(mode)
        staged.replace(path)
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def atomic_write_json(path: Path, value) -> None:
    """Serialize before opening a staged file; reject nonfinite JSON numbers."""
    text = json.dumps(value, indent=2, allow_nan=False) + "\n"
    with atomic_text_writer(path) as stream:
        stream.write(text)
