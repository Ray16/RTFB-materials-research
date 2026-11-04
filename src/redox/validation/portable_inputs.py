"""Read exact historical input bytes from an explicit portable mapping.

Historical bindings and currently executing producer files stay distinct.
Missing mappings never fall back to arbitrary original filesystem paths.
"""
import json
from pathlib import Path

from redox.core.hashing import sha256_file


class PortableInputs:
    def __init__(self, manifest_path, expected, analysis_sha256):
        self.manifest_path=Path(manifest_path).resolve()
        self.folder=self.manifest_path.parent
        manifest=json.loads(self.manifest_path.read_text())
        self.bindings=manifest['bindings'];self.local=set()
        if (manifest.get('schema')!=1 or manifest['analysis_sha256']!=analysis_sha256
                or {p:r['sha256'] for p,r in self.bindings.items()}!=expected):
            raise ValueError('Portable inputs must retain every exact historical binding')
        for row in self.bindings.values():
            name=row['file'];path=self.folder/name
            if (Path(name).name!=name or name in ('.','..') or path.is_symlink()
                    or not path.resolve().is_relative_to(self.folder)):
                raise ValueError('Unsafe portable input path')

    def allow_local(self, paths):
        """Authorize only explicit current input/producer paths supplied by caller."""
        self.local.update(str(Path(p).resolve()) for p in paths)

    def verify(self, original, expected):
        original=str(original);row=self.bindings.get(original)
        if row is not None and row['sha256']==expected:
            path=self.folder/row['file']
        elif str(Path(original).resolve()) in self.local:
            path=Path(original)
        else:
            raise ValueError('No exact portable input binding: '+original)
        if sha256_file(path)!=expected:
            raise ValueError('Portable input bytes changed: '+original)
        return path

    def path(self, original):
        row=self.bindings.get(str(original))
        if row is None:raise ValueError('Historical path is absent from portable inputs')
        return self.verify(original,row['sha256'])

    def digest(self, original):
        return sha256_file(self.path(original))
