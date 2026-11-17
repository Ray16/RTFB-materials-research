"""Read immutable evidence through one explicit checkout-prefix relocation.

Stored paths and hashes are never rewritten. An old checkout path is read only
from its declared archive location; a missing archive file cannot fall back to
an unrelated file now occupying the old checkout. External dependencies retain
their exact paths and must still pass their original content hashes.
"""
from pathlib import Path, PurePosixPath

from redox.core.hashing import sha256_file


class ArchivedSources:
    def __init__(self, prior_root, archive_root, *, relative_to_prior=False):
        prior = PurePosixPath(str(prior_root))
        if not prior.is_absolute() or '..' in prior.parts or prior == PurePosixPath('/'):
            raise ValueError('Declare an absolute original checkout prefix')
        self.prior = prior
        self.relative_to_prior = relative_to_prior
        self.archive = Path(archive_root).resolve(strict=True)
        if not self.archive.is_dir() or self.archive == Path(str(prior)).resolve():
            raise ValueError('A distinct existing archive directory is required')
        self._cache = {}
        self.receipts = {}

    def path(self, original):
        original = PurePosixPath(str(original))
        if '..' in original.parts:
            raise ValueError('Evidence paths must be absolute without parent traversal')
        if not original.is_absolute():
            if not self.relative_to_prior or original == PurePosixPath('.'):
                raise ValueError('Relative evidence requires an explicit original working-directory binding')
            original = self.prior/original
        try:
            relative = original.relative_to(self.prior)
        except ValueError:
            return Path(str(original))
        if not relative.parts:
            raise ValueError('Evidence must name a file below the original checkout')
        return self.archive.joinpath(*relative.parts)

    @staticmethod
    def _signature(path):
        value = path.stat()
        return (value.st_dev, value.st_ino, value.st_size,
                value.st_mtime_ns, value.st_ctime_ns)

    def verify(self, original, expected, *, force=False):
        if (not isinstance(expected, str) or len(expected) != 64
                or any(c not in '0123456789abcdef' for c in expected)):
            raise ValueError('A lowercase SHA-256 content digest is required')
        path = self.path(original)
        before = self._signature(path)
        cached = self._cache.get(str(path))
        actual = cached[1] if not force and cached and cached[0] == before else sha256_file(path)
        if self._signature(path) != before:
            raise ValueError('Evidence changed while hashing: '+str(path))
        if actual != expected:
            raise ValueError('Archived evidence content differs: '+str(original))
        self._cache[str(path)] = (before, actual)
        self.receipts[str(original)] = dict(relocated_path=str(path), sha256=actual)
        return path

    def verify_all(self, sources):
        if not isinstance(sources, dict) or not sources:
            raise ValueError('Nonempty source bindings required')
        for original, expected in sources.items():
            self.verify(original, expected)

    def recheck(self):
        for original, record in list(self.receipts.items()):
            self.verify(original, record['sha256'], force=True)
