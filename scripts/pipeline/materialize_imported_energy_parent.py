#!/usr/bin/env python
"""Make an identical local report alias for an imported energy receipt.

Native reports name the canonical local energy-report path. The imported
receipt still identifies its original source; this alias makes that path
resolvable without changing the source, receipt, plan or running calculation.
"""
import argparse
import json
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint


def materialize(campaign):
    records, invalid = api.admitted_records(campaign, "energy")
    if invalid or not records:
        raise ValueError("Checked imported electronic records are required")
    aliases = []
    for receipt in sorted((campaign / "energy/done").glob("*.json")):
        original = json.loads(receipt.read_text())
        if not original.get("reused_checked_energy"):
            raise ValueError("Only explicitly imported energy receipts are in scope")
        record = verify_state_checkpoint(receipt)
        source = Path(original["status_path"])
        destination = campaign / "energy/tasks" / record["task"]["key"] / "status.json"
        if source.resolve() == destination.resolve():
            continue
        if destination.exists():
            if destination.read_bytes() != source.read_bytes():
                raise ValueError("Preserve the conflicting canonical energy report")
        else:
            destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("xb") as output:
                output.write(source.read_bytes())
        if api.digest(destination) != original["status_sha256"]:
            raise ValueError("Imported parent alias differs from the original report")
        aliases.append(dict(source=str(source), alias=str(destination), sha256=api.digest(destination)))
    proof = dict(status="identical_imported_energy_parent_aliases_materialized", checked=time.time(),
        aliases=aliases, producer_sha256=api.digest(__file__), receipt_plan_and_source_changed=False,
        new_energy_calculations=0, new_native_G=0)
    atomic_write_json(campaign / "energy/imported_parent_aliases.json", proof)
    return proof


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign", type=Path)
    args = p.parse_args()
    print(json.dumps(materialize(args.campaign.resolve())))


if __name__ == "__main__":
    main()
