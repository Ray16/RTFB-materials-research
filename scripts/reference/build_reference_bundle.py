"""Freeze the pinned calibration evidence the pipeline needs into data/reference/.

The scorecard's uncertainties and the Fc/Fc+ reference are not constants: they are
re-derived live from raw records (DFT, active single points, native solution thermal,
spin evidence) of ferrocene and the grounded MeCN benchmark molecules, and the native
thermal gates also re-hash the certificate evidence cited in config/*.json. A fresh
checkout has none of these, so this script copies exactly that set from the
independently checked review's input inventory (review-time sha256 of every record,
recovering repaired files from their preserved superseded copies) into one deterministic archive plus a manifest
(per-file sha256, source, manifest rows, expected uncertainty evidence).

  python scripts/reference/build_reference_bundle.py --source ~/archive/0_redox \\
      --review results/reviews/2026-10-19/aryl_hydrogen_scope_tables

Seed a checkout with `redox seed-reference`; the bundle never recomputes anything.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAME = "validated_2026-10-19"
CALC_TREES = ("dft", "uma", "spincheck")
TABLES = ("orop_benchmark.csv", "d3tales_reorg_validation_prod/comparison.csv")
EVIDENCE_ROOTS = ("calcs/", "results/")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def certificate_evidence(source: Path) -> list[str]:
    """Every calcs/ or results/ file cited by a config certificate that exists in the source."""
    found = set()
    for config in sorted((ROOT / "config").glob("*.json")):
        for value in set(re.findall(r'"([^"\s]{4,400})"', config.read_text())):
            relative = value.split("/homes/rzhu/0_redox/", 1)[-1]
            if relative.startswith(EVIDENCE_ROOTS) and (source / relative).is_file():
                found.add(relative)
    return sorted(found)


def _owner(relative: str):
    parts = relative.split("/")
    if parts[0] == "library" and len(parts) > 2:
        return parts[1]
    if parts[0] == "calcs" and len(parts) > 3 and parts[1] in CALC_TREES:
        return parts[2]
    return None


def review_records(source: Path, ids: list[str], inventory: dict, source_root: str) -> dict[str, bytes]:
    """Exactly the record files the review read for `ids`, at their review-time sha256.

    Records later replaced by a repair are recovered from the preserved superseded copies
    (searched by content hash); a file that cannot be recovered aborts the build."""
    prefix = source_root.rstrip("/") + "/"
    wanted = {path[len(prefix):] if path.startswith(prefix) else path: digest
              for path, digest in inventory.items()}
    wanted = {rel: digest for rel, digest in wanted.items()
              if _owner(rel) in ids and not rel.endswith(".lock")}
    members, missing = {}, {}
    for relative, digest in sorted(wanted.items()):
        path = source / relative
        data = path.read_bytes() if path.is_file() else None
        if data is not None and sha256(data) == digest:
            members[relative] = data
        else:
            missing[digest] = relative
    if missing:
        roots = [source / "calcs" / tree / gid for tree in CALC_TREES for gid in {_owner(r) for r in missing.values()}]
        roots += [d for gid in {_owner(r) for r in missing.values()}
                  for d in (source / "results" / "validation").glob(f"*/*/{gid}")]
        for root in roots:
            for path in root.rglob("*") if root.is_dir() else ():
                if missing and path.is_file():
                    data = path.read_bytes()
                    relative = missing.pop(sha256(data), None)
                    if relative is not None:
                        members[relative] = data
        if missing:
            raise SystemExit(f"cannot recover review-time content for {sorted(missing.values())[:10]}")
    return members


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True, help="checked source tree (original checkout)")
    ap.add_argument("--review", required=True, help="source-relative review directory with the validated tables")
    ap.add_argument("--source-root", default="/nfs/lambda_stor_01/homes/rzhu/0_redox",
                    help="absolute checkout path the records were produced at")
    ap.add_argument("--out", type=Path, default=ROOT / "data" / "reference")
    args = ap.parse_args(argv)
    source, review = args.source.resolve(), args.source.resolve() / args.review

    expected = json.loads((review / "uncertainty_evidence.json").read_text())
    # E deg calibration events, disproportionation calibration molecules, and the Fc reference.
    stability = {row["id"] for row in csv.DictReader(open(review / "stability_validation.csv", newline=""))}
    ids = sorted({event["id"] for event in expected["grounded_redox_events"]} | stability | {"ferrocene"})
    rows = [row for row in csv.DictReader(open(source / "library" / "manifest.csv", newline=""))
            if row["id"] in ids]
    if {row["id"] for row in rows} != set(ids):
        raise SystemExit(f"manifest rows missing for {sorted(set(ids) - {r['id'] for r in rows})}")

    inventory = json.loads((review / "input_inventory.json").read_text())
    members = review_records(source, ids, inventory, args.source_root)
    for relative in certificate_evidence(source):
        members.setdefault(relative, (source / relative).read_bytes())
    for table in TABLES:
        members[f"results/{table}"] = (review / table).read_bytes()

    args.out.mkdir(parents=True, exist_ok=True)
    archive = args.out / f"{NAME}.tar.xz"
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:xz", preset=6) as tar:
        for relative in sorted(members):
            info = tarfile.TarInfo(relative)
            info.size, info.mtime, info.mode = len(members[relative]), 0, 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(members[relative]))
    archive.write_bytes(buffer.getvalue())

    manifest = dict(
        schema=1, name=NAME,
        source_tree=str(source), source_root=args.source_root,
        source_git_commit=_git_head(source), review=args.review,
        purpose="pinned calibration evidence: Fc/Fc+ reference, grounded MeCN E deg and disproportionation benchmark records",
        archive=archive.name, archive_sha256=sha256(archive.read_bytes()),
        ids=ids, manifest_rows=rows,
        files={relative: sha256(data) for relative, data in sorted(members.items())},
        expected_uncertainty_evidence=expected,
        expected_fc_abs_V=json.loads((ROOT / "config" / "project.json").read_text())["reference"]["fc_abs_computed_V"],
    )
    (args.out / f"{NAME}.manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    print(f"{archive}: {len(members)} files, {archive.stat().st_size / 1e6:.1f} MB; ids={len(ids)}")


def _git_head(source: Path) -> str | None:
    head = source / ".git" / "HEAD"
    try:
        ref = head.read_text().strip()
        if ref.startswith("ref: "):
            return (source / ".git" / ref[5:]).read_text().strip()
        return ref
    except OSError:
        return None


if __name__ == "__main__":
    sys.exit(main())
