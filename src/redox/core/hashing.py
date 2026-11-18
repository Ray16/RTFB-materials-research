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


def verify_sources(mapping):
    """Raise if any path's current bytes no longer match its recorded digest."""
    for path, expected in mapping.items():
        if sha256_file(path) != expected:
            raise RuntimeError(f"Source or input changed: {path}")


_FILE_HASH_CACHE: dict = {}


def sha256_file_cached(path):
    """`sha256_file`, memoized per (path, mtime, size) for this process.

    Provenance checks re-hash the same immutable, multi-GiB QM artifacts (Hessian
    .npy, SCF.h5, response blocks) many times while one stage reads properties for
    several states. Hashing is bounded-memory (streamed) and, once computed, reused
    until the file changes on disk. Missing files raise, exactly as reading would.
    """
    p = Path(path)
    st = p.stat()
    key = (str(p), st.st_mtime_ns, st.st_size)
    value = _FILE_HASH_CACHE.get(key)
    if value is None:
        value = sha256_file(p)
        _FILE_HASH_CACHE[key] = value
    return value
