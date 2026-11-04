"""Relocate mutable result pointers without rewriting scientific evidence."""
from pathlib import Path


def relocate_index(index, relocation, root):
    """Return a copied index and a list of verified pointer changes.

    Only exact old directory prefixes from the explicit relocation record may
    change. Every new target must exist inside the canonical project. Historical
    evidence files are never opened or rewritten. Live/deferred roots have no
    mapping and remain untouched. Existing old paths must be the recorded links.
    """
    root = Path(root).resolve()
    mappings = []
    for row in relocation['moved']:
        old, new = Path(row['old']), Path(row['new'])
        if not old.is_absolute() or not new.is_absolute() or '..' in old.parts or '..' in new.parts:
            raise ValueError('Relocation prefixes must be absolute without parent traversal')
        if not new.is_dir() or not new.resolve().is_relative_to(root):
            raise ValueError('Relocation destination must exist inside the canonical project')
        if old.exists() or old.is_symlink():
            if not old.is_symlink() or old.resolve() != new.resolve():
                raise ValueError('An old location is not its recorded compatibility link')
        if any(old == prior or old.is_relative_to(prior) or prior.is_relative_to(old)
               for prior, _ in mappings):
            raise ValueError('Relocation prefixes must not overlap')
        mappings.append((old, new))
    changes = []

    def visit(value, location):
        if isinstance(value, dict):
            return {key: visit(item, location + [key]) for key, item in value.items()}
        if isinstance(value, list):
            return [visit(item, location + [i]) for i, item in enumerate(value)]
        if isinstance(value, str):
            for old, new in mappings:
                if value == str(old) or value.startswith(str(old) + '/'):
                    relative = Path(value).relative_to(old)
                    target = new / relative
                    if '..' in relative.parts or not target.exists() or not target.resolve().is_relative_to(root):
                        raise ValueError('Canonical pointer target is missing or escapes the project: ' + value)
                    changes.append(dict(field=location, original=value, canonical=str(target)))
                    return str(target)
        return value

    return visit(index, []), changes


def relocated_input_path(path: str | Path, relocation: dict, root: str | Path) -> Path:
    """Resolve one input through an explicit, checked workspace mapping.

    Unmapped inputs keep their declared location. This changes no source record;
    callers must still verify the input's expected content hash before use.
    """
    resolved, _ = relocate_index(str(path), relocation, root)
    return Path(resolved)
