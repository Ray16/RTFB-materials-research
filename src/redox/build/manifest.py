"""One writer contract for candidate, standalone and validation libraries."""
from contextlib import contextmanager
import csv
import fcntl
from pathlib import Path

from redox.core.io import atomic_text_writer

COLS = ["id", "name", "family", "state", "charge", "mult_hint",
        "n_e", "n_conf", "solv_preopt", "smiles"]


@contextmanager
def locked_manifest(path: Path):
    """Hold across reading, building structures and replacing the merged manifest."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with (path.parent / ".manifest.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def read_rows(path: Path) -> list[dict]:
    path = Path(path)
    if not path.exists():
        return []
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows) -> None:
    """Replace a manifest atomically. The caller must hold ``locked_manifest``."""
    with atomic_text_writer(path, newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLS)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in COLS} for row in rows)
