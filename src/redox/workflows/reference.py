"""Seed a checkout with the pinned calibration evidence in data/reference/.

Every Fc-referenced potential and every scorecard uncertainty is re-derived live from raw
records (ferrocene; grounded MeCN benchmark molecules) and from certificate evidence that
the native thermal gates re-hash. This installs exactly those immutable files into
library/, calcs/ and results/, verifying each sha256, and never overwrites a differing
file. `--verify` then recomputes the Fc reference and the uncertainty evidence through the
production gates and requires an exact match with the validated values.

  python -m redox seed-reference            # install (idempotent) and verify
  python -m redox seed-reference --check    # report only, write nothing
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import tarfile
from pathlib import Path

from redox.core.common import ROOT

REFERENCE = ROOT / "data" / "reference"
DEFAULT = "validated_2026-10-19"


def load(name: str = DEFAULT, reference: Path = REFERENCE):
    manifest = json.loads((reference / f"{name}.manifest.json").read_text())
    data = (reference / manifest["archive"]).read_bytes()
    if hashlib.sha256(data).hexdigest() != manifest["archive_sha256"]:
        raise RuntimeError(f"{manifest['archive']}: archive sha256 differs from its manifest")
    members = {}
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:xz") as tar:
        for info in tar.getmembers():
            if not info.isfile() or info.name.startswith("/") or ".." in Path(info.name).parts:
                raise RuntimeError(f"unsafe archive member {info.name!r}")
            members[info.name] = tar.extractfile(info).read()
    if set(members) != set(manifest["files"]):
        raise RuntimeError("archive members differ from the manifest file list")
    for name, content in members.items():
        if hashlib.sha256(content).hexdigest() != manifest["files"][name]:
            raise RuntimeError(f"{name}: sha256 differs from the manifest")
    return manifest, members


def plan(manifest, members, root: Path = ROOT):
    """Classify every file as present (identical), missing, or conflicting."""
    present, missing, conflicts = [], [], []
    for name, content in sorted(members.items()):
        target = root / name
        if not target.exists():
            missing.append(name)
        elif target.read_bytes() == content:
            present.append(name)
        else:
            conflicts.append(name)
    return present, missing, conflicts


def merge_manifest_rows(rows, root: Path = ROOT):
    """Add pinned rows; an existing row for the same (id, state) must match exactly."""
    from redox.structures.manifest import COLS, locked_manifest, read_rows, write_rows
    path = root / "library" / "manifest.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with locked_manifest(path):
        existing = read_rows(path) if path.exists() else []
        index = {(r["id"], r["state"]): r for r in existing}
        added = []
        for row in rows:
            row = {key: str(row.get(key, "")) for key in COLS}
            current = index.get((row["id"], row["state"]))
            if current is None:
                added.append(row)
            elif {k: str(current.get(k, "")) for k in COLS} != row:
                raise RuntimeError(f"manifest row {row['id']}/{row['state']} differs from the pinned reference")
        if added:
            write_rows(path, existing + added)
    return len(added)


def live_evidence(manifest):
    """Recompute Fc and the uncertainty evidence from the installed records."""
    from redox.core.common import free_energy, read_result
    from redox.properties.potentials import fc_live
    from redox.screening.scorecard import derive_uncertainties
    reference = fc_live()
    rows = []
    for event in manifest["expected_uncertainty_evidence"]["grounded_redox_events"]:
        ox, red = event["event"].split("->")
        (g_ox, q_ox), (g_red, q_red) = (free_energy(read_result(event["id"], s)) for s in (ox, red))
        value = (round(g_ox - g_red - reference, 3)
                 if None not in (g_ox, g_red, reference) and q_ox == q_red == "ok" else None)
        rows.append(dict(id=event["id"], event=event["event"], E_vs_Fc_V=value))
    return reference, derive_uncertainties(rows)


# Derived by the finalization tables (redox.properties.stability -> redox.quality.stability)
# from these same records, so it is checked against the validated value after finalization.
TABLE_DERIVED = ("sigma_disp",)


def verify(manifest):
    reference, evidence = live_evidence(manifest)
    expected = {k: v for k, v in manifest["expected_uncertainty_evidence"].items() if k not in TABLE_DERIVED}
    problems = []
    if reference is None or abs(reference - manifest["expected_fc_abs_V"]) > 1e-9:
        problems.append(f"Fc reference {reference!r} != pinned {manifest['expected_fc_abs_V']!r}")
    for key, value in expected.items():
        if evidence.get(key) != value:
            problems.append(f"{key}: live {str(evidence.get(key))[:120]} != validated {str(value)[:120]}")
    return problems


def compare_results(manifest, results: Path = ROOT / "results"):
    """Finalized sigmas derive only from the pinned benchmark records, so registering new
    candidates must leave them unchanged; any difference is a calibration change to review."""
    path = results / "uncertainty_evidence.json"
    if not path.is_file():
        print(f"[reference] {path} missing: run the finalization tables first")
        return 2
    current = json.loads(path.read_text())
    problems = [f"{key}: finalized {str(current.get(key))[:120]} != validated {str(value)[:120]}"
                for key, value in manifest["expected_uncertainty_evidence"].items() if current.get(key) != value]
    if problems:
        print("[reference] CALIBRATION CHANGED vs pinned evidence:\n  " + "\n  ".join(problems))
        return 2
    print(f"[reference] finalized uncertainty evidence equals the validated {manifest['name']} evidence "
          f"({len(manifest['expected_uncertainty_evidence'])} fields)")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default=DEFAULT)
    ap.add_argument("--check", action="store_true", help="report only; write nothing")
    ap.add_argument("--no-verify", action="store_true", help="skip the live recomputation check")
    ap.add_argument("--compare-results", action="store_true",
                    help="after finalization: results/uncertainty_evidence.json must equal the validated evidence")
    args = ap.parse_args(argv)
    manifest, members = load(args.name)
    if args.compare_results:
        return compare_results(manifest)
    present, missing, conflicts = plan(manifest, members)
    print(f"[reference] {args.name}: {len(present)} present, {len(missing)} missing, {len(conflicts)} conflicting")
    if conflicts:
        for name in conflicts[:20]:
            print(f"  CONFLICT {name}")
        raise SystemExit("Refusing to overwrite differing files; resolve them before seeding")
    if args.check:
        return 1 if missing else 0
    for name in missing:
        target = ROOT / name
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(target.name + ".seeding")
        temporary.write_bytes(members[name])
        temporary.replace(target)
    added = merge_manifest_rows(manifest["manifest_rows"])
    print(f"[reference] installed {len(missing)} files; {added} manifest rows added")
    if str(ROOT) != manifest["source_root"]:
        print(f"[reference] NOTE: checkout {ROOT} differs from the records' origin {manifest['source_root']}; "
              "native thermal certificates cite that origin, so admission needs "
              "redox.quality.relocated_native.relocated_native_reader")
    if args.no_verify:
        return 0
    problems = verify(manifest)
    if problems:
        print("[reference] VERIFY FAILED:\n  " + "\n  ".join(problems))
        return 2
    print(f"[reference] verified: Fc {manifest['expected_fc_abs_V']} V and "
          f"{len(manifest['expected_uncertainty_evidence']) - len(TABLE_DERIVED)} record-derived uncertainty "
          f"fields reproduce exactly ({', '.join(TABLE_DERIVED)} is checked after finalization)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
