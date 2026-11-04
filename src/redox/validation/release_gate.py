"""Admission checks for an archive whose streamed byte verification has passed."""
from pathlib import Path


def verified_archive_entry(verification, raw_manifest, manifest_sha256, inventory, bundle):
    """Validate a trusted verifier receipt against its inventory and local file.

    The verifier has already hashed the full archive and every member. This gate
    checks that its receipt belongs to the exact manifest/inventory being added;
    it does not repeat a large archive read or establish physical qualification.
    """
    objects = raw_manifest['objects']
    digest = verification.get('archive_sha256')
    if (verification.get('status') != 'all_streamed_archive_members_verified'
            or verification.get('manifest_sha256') != manifest_sha256
            or raw_manifest.get('inventory_sha256') != inventory['inventory_sha256']
            or verification.get('objects') != len(objects)
            or len(objects) != inventory['archive_objects']
            or sum(objects.values()) != inventory['research_evidence_bytes']
            or verification.get('members_verified') != len(objects) + len(raw_manifest['support_files_sha256']) + 1
            or not isinstance(digest, str) or len(digest) != 64
            or any(c not in '0123456789abcdef' for c in digest)
            or verification.get('numerical_replay_performed') is not False):
        raise ValueError('Raw archive has no matching complete byte-verification receipt')
    archive = Path(verification['archive'])
    bundle = Path(bundle).resolve()
    if (archive.is_symlink() or not archive.is_file()
            or not archive.resolve().is_relative_to(bundle)
            or archive.stat().st_size != verification['archive_bytes']):
        raise ValueError('Verified archive is missing, outside the package or has changed size')
    return dict(id='full_declared_raw_evidence',file=str(archive.resolve().relative_to(bundle)),
                sha256=digest,bytes=verification['archive_bytes'],objects=len(objects),
                scope='Declared research-evidence bytes and mappings; inclusion adds no physical or experimental qualification')
