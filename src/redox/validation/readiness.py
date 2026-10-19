"""Manifest-based readiness, including states that have never produced any output."""
import argparse
import csv
import json
from pathlib import Path

from redox.core.common import ROOT, DFT, read_manifest, free_energy, apply_thermal_policy
from redox.build.geometry import smiles_to_xyz_map, geometry_stereo_status
from redox.core.protocol import load_record


def inspect_states(rows=None, root=DFT):
    rows = read_manifest() if rows is None else rows
    out = []
    for row in rows:
        gid, st = row["id"], row["state"]
        sd = Path(root) / gid / st
        issues = []
        result = sd / "result.json"
        rec = None
        if not result.exists() or not (sd / "opt.xyz").exists():
            issues.append("missing_geometry_or_result")
        else:
            try:
                r = json.loads(result.read_text())
                q, m = int(r["charge"]), int(r["mult"])
                if q != int(row["charge"]):
                    issues.append("charge_mismatch")
                if row.get("smiles"):
                    from rdkit import Chem
                    declared = Chem.AddHs(Chem.MolFromSmiles(row["smiles"]))
                    actual = Chem.MolFromXYZFile(str(sd / "opt.xyz"))
                    if actual is None or sorted(a.GetAtomicNum() for a in actual.GetAtoms()) != sorted(
                            a.GetAtomicNum() for a in declared.GetAtoms()):
                        issues.append("atom_composition_mismatch")
                    elif smiles_to_xyz_map(row["smiles"], sd / "opt.xyz") is None:
                        issues.append("declared_connectivity_mismatch")
                    elif geometry_stereo_status(row['smiles'], sd / 'opt.xyz') not in ('ok', 'not_applicable'):
                        issues.append('declared_stereochemistry_mismatch')
                    if actual is not None:
                        electrons = sum(a.GetAtomicNum() for a in actual.GetAtoms()) - q
                        if m < 1 or electrons < m - 1 or (electrons - (m - 1)) % 2:
                            issues.append("spin_parity_mismatch")
                from redox.qm.dft import _opt_protocol_mismatch
                if _opt_protocol_mismatch(r, q):
                    issues.append("stale_geometry_protocol")
                rec = load_record(sd, q, m)
                if rec is None:
                    issues.append("missing_or_invalid_active_sp")
                G, tq = free_energy(apply_thermal_policy(dict(r, e_smd_eV=0.0), sd))
                if G is None:
                    issues.append(tq if tq in ("missing_validated_native_thermal", "invalid_thermal_policy",
                                              "wrong_spin", "stale_spin_decision")
                                  else "missing_thermal")
                elif tq == "imag_modes":
                    issues.append("imaginary_thermal_modes")
            except (ValueError, KeyError, TypeError):
                issues.append("malformed_result")
        blocking = [v for v in issues if v != "imaginary_thermal_modes"]
        invalid = {"atom_composition_mismatch", "declared_connectivity_mismatch", "declared_stereochemistry_mismatch", "spin_parity_mismatch", "charge_mismatch", "wrong_spin"}
        out.append(dict(id=gid, state=st, charge=row["charge"],
                        status="INVALID" if invalid.intersection(issues) else (
                            "INCOMPLETE" if blocking else ("REVIEW" if issues else "ready")),
                        issues=";".join(issues)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", help="comma-separated molecule IDs")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--strict", action="store_true", help="exit nonzero unless every selected state is ready")
    ap.add_argument("--include-lineage", action="store_true", help="include hash-verified original repair geometries in the CSV; historical records stay unchanged")
    a = ap.parse_args()
    rows = read_manifest()
    if a.only:
        ids = set(a.only.split(","))
        unknown = ids - {r["id"] for r in rows}
        if unknown:
            ap.error(f"unknown IDs: {sorted(unknown)}")
        rows = [r for r in rows if r["id"] in ids]
    report = inspect_states(rows)
    fields=["id", "state", "charge", "status", "issues"]
    if a.include_lineage:
        from redox.validation.provenance import repair_source
        fields += ['repair_source_status','repair_source_geometry','repair_source_sha256']
        for item in report:
            try:
                raw=json.loads((DFT/item['id']/item['state']/'result.json').read_text())
            except (OSError,ValueError):
                raw={}
            item.update(repair_source(raw,ROOT))
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        with a.out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader(); w.writerows(report)
    for r in report:
        if r["status"] != "ready":
            print(f"{r['id']}/{r['state']}: {r['status']} {r['issues']}")
    print(f"{len(report)} declared states; " + ", ".join(
        f"{sum(r['status'] == s for r in report)} {s}" for s in ("ready", "REVIEW", "INCOMPLETE", "INVALID")))
    if a.include_lineage:
        print(f"Repair sources: {sum(r['repair_source_status']=='verified' for r in report)} direct, "
              f"{sum(r['repair_source_status']=='resolved_archived_source' for r in report)} verified archives, "
              f"{sum(r['repair_source_status']=='unverified' for r in report)} unresolved")
    if a.strict and any(r["status"] != "ready" for r in report):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
