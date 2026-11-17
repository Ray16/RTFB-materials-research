"""Bounded-memory file hashing for large calculation inventories."""
from pathlib import Path
import hashlib


def sha256_file(path, *, chunk_size=1024 * 1024):
    """Hash file bytes without loading a Hessian or checkpoint into memory."""
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()
