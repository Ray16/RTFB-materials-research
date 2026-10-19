"""Existing file dependencies named in immutable validation evidence."""
from pathlib import Path
import hashlib


def calculation_gid(path, root):
    """Identify molecule-owned calculation inputs; unknown/global paths stay global."""
    try:
        parts = Path(path).resolve().relative_to(Path(root).resolve() / 'calcs').parts
    except ValueError:
        return None
    return parts[1] if len(parts) >= 3 and parts[0] in ('dft', 'uma', 'spincheck') else None


def scoped_spin_summary(summary, gids, root):
    """Limit publication freshness to displayed molecules and their frozen evidence.

    Older snapshots lack per-molecule inventories and keep their broader checks.
    Scientific code, configuration and unclassified dependencies are always retained.
    """
    inventories = summary.get('spin_input_paths_by_gid')
    if inventories is None:
        return summary
    retained = set(gids) & set(summary['spin_input_gids'])
    result = dict(summary)
    result['spin_input_gids'] = sorted(retained)
    result['spin_input_paths'] = sorted({path for gid in retained for path in inventories[gid]})
    result['admitted_input_sha256'] = {
        path: digest for path, digest in summary.get('admitted_input_sha256', {}).items()
        if calculation_gid(path, root) in retained or calculation_gid(path, root) is None
    }
    return result


def input_changes_quarantined(summary):
    """Changed inputs may publish only after every dependent molecule is excluded."""
    changes = summary.get('inputs_changed_during_snapshot') or []
    if not changes:
        return True
    if summary.get('global_input_changed') is not False:
        return False
    retained = set(summary.get('spin_input_gids', []))
    scopes = summary.get('input_change_scope') or {}
    for path in changes:
        owners = set(scopes.get(path, []))
        if not owners or 'ferrocene' in owners or owners & retained:
            return False
    return True


def repair_source(record, root):
    """Resolve a recorded original geometry by hash, preserving historical bytes.

    Older promotion records retained a pointer to their replaced canonical
    geometry. Only an explicit archive whose bytes match the original hash can
    resolve that pointer; neither directory age nor a newer geometry is evidence.
    """
    if not isinstance(record, dict):
        return dict(repair_source_status='unverified', repair_source_geometry=None,
                    repair_source_sha256=None)
    source = record.get('minimum_repair_source')
    expected = record.get('minimum_repair_source_sha256')
    result = dict(repair_source_status='not_declared', repair_source_geometry=None,
                  repair_source_sha256=expected)
    if not source and not expected:
        return result
    result['repair_source_status'] = 'unverified'
    if not isinstance(expected, str) or len(expected) != 64 or any(c not in '0123456789abcdef' for c in expected):
        return result

    def matches(path):
        try:
            return hashlib.sha256(path.read_bytes()).hexdigest() == expected
        except OSError:
            return False

    if isinstance(source, str) and matches(Path(root) / source):
        result.update(repair_source_status='verified', repair_source_geometry=str(Path(root) / source))
        return result
    archive = record.get('legacy_thermal_archive')
    if isinstance(archive, str):
        candidate = (Path(root) / archive).parent / 'opt.xyz'
        if matches(candidate):
            result.update(repair_source_status='resolved_archived_source', repair_source_geometry=str(candidate))
    return result


def referenced_files(value, root):
    paths = set()
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, str) and (key == 'path' or key.endswith('_path')):
                candidate = Path(root) / item
                if candidate.is_file(): paths.add(candidate)
            if isinstance(item, dict) and key.endswith('sha256'):
                for name, digest in item.items():
                    if isinstance(digest, str) and len(digest) == 64 and '/' in name:
                        candidate = Path(root) / name
                        if candidate.is_file(): paths.add(candidate)
            paths.update(referenced_files(item, root))
    elif isinstance(value, list):
        for item in value: paths.update(referenced_files(item, root))
    return paths
